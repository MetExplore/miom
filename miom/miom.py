import numpy as np
import warnings
import mip
from abc import ABC, abstractmethod
from functools import wraps
from collections.abc import Iterable
from collections import defaultdict
from miom.mio import load_gem
from typing import NamedTuple
from enum import Enum, auto
from time import perf_counter

try:
    import picos as pc
    _picos_backend_available = True
except ImportError:
    _picos_backend_available = False


class Status(str, Enum):
    """Subset of PICOS solver status codes.

    This subset of states offers compatibility between Python-MIP and PICOS.
    Status codes can be obtained after solving an optimization problem with
    the property `status`. For example:

    ```python
    >>> import miom
    >>> miom.load('@iMM1865').steady_state().set_rxn_objective('BIOMASS_reaction').solve().status
    {'status': 'optimal', 'objective_value': 798.811, 'elapsed_seconds': 0.479}
    ```
    """
    OPTIMAL = "optimal",
    FEASIBLE = "feasible"
    INFEASIBLE = "infeasible",
    UNBOUNDED = "unbounded",
    ILLPOSED = "illposed",
    EMPTY = "empty",
    PREMATURE = "premature",
    UNKNOWN = "unknown"


class Solvers(str, Enum):
    """Solvers supported by the miom module.
    
    Please refer to https://picos-api.gitlab.io/picos/introduction.html to see
    the list of supported solvers using the PICOS backend. The Python-MIP backend
    only supports the GUROBI and CBC solvers.

    Note that in some cases, the Python-MIP backend (GUROBI/CBC) might be faster
    setting up the problem than the PICOS backend since it has less overhead.

    Attributes:
        GUROBI_PYMIP (str): Recommended solver for most problems (LP, MIP).
            It uses the Python-MIP backend. Note that GUROBI is a commercial
            software which require a license. Free academic licenses are 
            also available at https://gurobi.com/free/.
        GUROBI (str): For using GUROBI with the PICOS backend instead of Python-MIP.
        COIN_OR_CBC (str): Free LP/MIP solver with good performance, provided with Python-MIP.
        CPLEX (str): Commercial solver with a performance similar to GUROBI. Only
            supported by the PICOS backend.
        SCIP (str): Free academic solver, supported by the PICOS backend.
        GLPK (str): Open source solver, supported by the PICOS backend.
        MOSEK (str): Commercial solver, supported by the PICOS backend.
    """
    GUROBI_PYMIP = "gurobi_pymip",
    GUROBI = "gurobi",
    COIN_OR_CBC = "cbc",
    CPLEX = "cplex",
    GLPK = "glpk",
    SCIP = "scip",
    CVXOPT = "cvxopt",
    MOSEK = "mosek"


class _ReactionType(Enum):
    RH_POS = auto(),
    RH_NEG = auto(),
    RL = auto()


class ExtractionMode(str, Enum):
    """Method to select the subnetwork to be extracted.

    This mode works with the method [select_subnetwork][miom.miom.BaseModel] 
    only after a subset selection problem was succesfully solved.
    If the model is configured as a subset selection problem, the user can
    extract a subnetwork after the method solve() was called. The selection
    of the subnetwork can be done using the indicator variables or the flux
    variables.

    For more information, please read the documentation of the method
    [subset_selection][miom.miom.BaseModel]

    Attributes:
        ABSOLUTE_FLUX_VALUE (str): Use a decision criterion based on the
            value of the fluxes. For example, by selecting all the reactions
            that have an absolute flux value above a certain threshold.

        INDICATOR_VALUE (str): Use the binary indicator variables to select
            the subnetwork. Binary indicator variables are created for a
            subset selection problem (after calling [subset_selection][miom.miom.BaseModel]).
            Two indicators are created for each positive weighted and reversible reaction,
            to indicate if there is a non-zero positive flux or negative flux respectively.
            A single binary indicator variable is created for each negative weighted reaction.
            In this case, an indicator value of 1 indicates that the reaction was succesfully
            removed, and 0 otherwise. You can use the value of the indicator variables to
            select a subnetwork after solving. For example, if all the reactions have a
            negative weight (since the goal is, for example, to minimize the number of 
            reactions subject to some other constraints) and you want to select only the
            reactions that were not removed after solving, you can use the indicator
            variables with a value of 0.

            Usage example:
            ```python
                m.
                steady_state().
                add_constraints(...).
                # Convert to a subset selection, where each reaction
                # has a weight of -1.
                subset_selection(-1).
                # Solve the optimization problem
                .solve()
                # Get the subnetwork selecting the reactions
                # with an indicator value below 0.5. Since variables
                # are binary, it basically selects all the reactions
                # associated with a binary indicator value of 0.
                # This corresponds to the reactions that were not
                # removed after solving the problem (since all the
                # reactions have negative weights, and their binary
                # indicator variables are 1 if they were successfully
                # removed, and 0 if they could not be removed).
                .select_subnetwork(
                    mode=ExtractionMode.INDICATOR_VALUE,
                    comparator=Comparator.LESS_OR_EQUAL,
                    value=0.5
                ).network
            ```
    """
    ABSOLUTE_FLUX_VALUE = "flux_value",
    INDICATOR_VALUE = "indicator_value",
    REACTION_ACTIVITY = "reaction_activity"


class Comparator(str, Enum):
    """Comparator enum to use with [select_subnetwork()][miom.miom.BaseModel]

    Attributes:
        GREATER_OR_EQUAL (str): Select those variables 
            with a value greater or equal than the one provided.
        LESS_OR_EQUAL (str): Select those variables 
            with a value less or equal than the one provided.
    """
    GREATER_OR_EQUAL = "geq",
    LESS_OR_EQUAL = "leq"


_RxnVar = NamedTuple(
    "RxnVar", [
        ("index", int),
        ("id", int),
        ("lb", float),
        ("ub", float),
        ("cost", float),
        ("type", _ReactionType)
    ])

def _get_reversible_vars(weighted_rxns):
    d = defaultdict(list)
    for i, rxn in enumerate(weighted_rxns):
        if rxn.type == _ReactionType.RH_POS or rxn.type == _ReactionType.RH_NEG:
            d[rxn.index].append(i)
    return [l for l in d.values() if len(l) > 1]


def _get_rxn_var_data(weighted_rxns, rxn_type):
    idx = [
        (i, rd.index, rd.lb, rd.ub) 
        for i, rd in enumerate(weighted_rxns) 
        if rd.type is rxn_type
    ]
    if len(idx) > 0:
        return [np.array(x) for x in list(zip(*idx))]
    return None


def load(network, solver=None):
    """
    Create a MIOM optimization model for a given solver.
    If the solver is Coin-OR CBC, an instance of PythonMipModel is used
    (which uses the Python-MIP lib as the backend). Otherwise, a PicosModel
    is created (which uses PICOS as the backend).

    Example:
        Example of how to perform FBA to maximize flux through the
        `BIOMASS_reaction` in the iMM1865 model:

        ```python
        >>> import miom
        >>> network = miom.load_gem("@iMM1865")
        >>> flx = (miom
                    .load(network)
                    .steady_state()
                    .set_rxn_objective("BIOMASS_reaction")
                    .solve(verbosity=1)
                    .get_fluxes('BIOMASS_reaction'))
        ```

    Args:
        network (miom_network or str): A miom metabolic network or a valid file that
            can be imported with [load_gem()][miom.mio.load_gem].
        solver (Solver, optional): The solver to be used. Defaults to Solver.GLPK.

    Returns:
        BaseModel: A BaseModel object, which can be PythonMipModel if CBC solver is
            used, or a PicosModel otherwise.
    """
    if solver is None:
        solver = Solvers.COIN_OR_CBC
        if _picos_backend_available:
            solvers = pc.solvers.available_solvers()
            if "gurobi" in solvers:
                solver = Solvers.GUROBI
            elif "cplex" in solvers:
                solver = Solvers.CPLEX

    solver = str(solver.value) if isinstance(solver, Enum) else str(solver)
    if isinstance(network, str):
        network = load_gem(network)
    if solver == 'cbc':
        return PythonMipModel(miom_network=network, solver_name=solver)
    if solver == 'gurobi_pymip':
        return PythonMipModel(miom_network=network, solver_name='gurobi')
    if _picos_backend_available:
        return PicosModel(miom_network=network, solver_name=solver)
    else:
        raise Exception("""PICOS is not installed. Please install it with pip, 
                        or reinstall miom with the right options, for example with:
                        'pip install miom[all]'.""")


class _Variables(ABC):
    def __init__(self, flux_vars=None, indicator_vars=None, assigned_reactions=None):
        self._flux_vars = flux_vars
        self._indicator_vars = indicator_vars
        self._assigned_reactions = assigned_reactions

    @property
    def indicators(self):
        return self._indicator_vars

    @property
    def fluxvars(self):
        return self._flux_vars

    @property
    @abstractmethod
    def flux_values(self):
        pass

    @property
    @abstractmethod
    def indicator_values(self):
        pass

    @property
    def reaction_activity(self):
        """Returns a list indicating whether a reaction is active or not.

        It uses the value of the indicator variables of the subset selection
        problem to indicate whether a reaction is active (has positive or
        negative flux, in which case the value is 1 or -1. A value of None 
        indicates that the reaction has no associated indicator variable. 
        A value of `nan` indicates that the value of the reaction is 
        inconsistent after solving the subset selection problem.

        Returns:
            list: list of the same length as the number of reactions in the model.
        """
        activity = [None] * len(self.fluxvars)
        values = self.indicator_values
        if values is None:
            raise ValueError("The indicator values are not set. This means that \
                the problem is not a MIP (subset_selection method was not called). \
                If you want to select a subset of reactions based on a flux value, \
                use the method select_subnetwork or obtain_subnetwork instead.")
        for i, rxn in enumerate(self.assigned_reactions):
            curr_value = activity[rxn.index]
            if rxn.type == _ReactionType.RH_POS or rxn.type == _ReactionType.RH_NEG:
                v = 0 if curr_value is None else curr_value
                v += values[i]
                activity[rxn.index] = v
            elif rxn.type == _ReactionType.RL:
                if curr_value is not None:
                    raise ValueError("Multiple indicator variables for the same RL type reaction")
                activity[rxn.index] = 1 - values[i]
        # Replace inconsistent values (can happen due to numerical issues)
        for i in range(len(activity)):
            if activity[i] is not None and activity[i] > 1:
                activity[i] = float('nan')
        # Add sign
        for i, rxn in enumerate(self.assigned_reactions):
            if rxn.type == _ReactionType.RH_NEG and values[i] > 0:
                activity[rxn.index] = -1 * activity[rxn.index]
        return activity


    @property
    def assigned_reactions(self):
        return self._assigned_reactions

    def values(self):
        return self.flux_values, self.indicator_values


class _PicosVariables(_Variables):
    def __init__(self):
        super().__init__()

    @property
    def flux_values(self):
        if self.fluxvars is None:
            return None
        arr = np.array(self.fluxvars.value).squeeze()
        if len(arr) == 1:
            return arr[0]
        return arr

    @property
    def indicator_values(self):
        if self.indicators is None:
            return None
        arr = np.array(self.indicators.value).squeeze()
        if len(arr) == 1:
            return arr[0]
        return arr


class _PythonMipVariables(_Variables):
    def __init__(self):
        super().__init__()

    @property
    def flux_values(self):
        if self.fluxvars is not None:
            return np.array([v.x for v in self.fluxvars])
        return None

    @property
    def indicator_values(self):
        if self.indicators is not None:
            return np.array([v.x for v in self.indicators])
        return None


def _composable(fn):
    """Annotation for methods that return the instance itself to enable chaining.

    If a method `my_method` is annotated with @_composable, a method called `_my_method`
    is expected to be provided by a subclass. Parent method `my_method` is called first
    and the result is passed to the child method `_my_method`.

    """
    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        # Invoke first the original method
        result = fn(self, *args, **kwargs)
        # If the result is a dict, update the kwargs
        if isinstance(result, dict):
            kwargs.update(result)
        else:
            kwargs['_parent_result'] = result
        if result is False:
            return self
        # Find subclass implementation
        fname = '_' + fn.__name__
        if not hasattr(self, fname):
            raise ValueError(f'Method "{fn.__name__}()" is marked as @_composable '
                             f'but the expected implementation "{fname}()" was not provided '
                             f'by {self.__class__.__name__}')
        func = getattr(self, fname)
        result = func(*args, **kwargs)
        if isinstance(result, BaseModel):
            return result
        return self

    return wrapper


def _weighted_rxns(R, weights=None):
    rxn_data = []
    for i, rxn in enumerate(R):
        id, lb, ub = rxn['id'], rxn['lb'], rxn['ub']
        w = weights[i] if weights is not None else 0
        if w > 0 and ub > 0:
            rxn_data.append(_RxnVar(i, id, lb, ub, w, _ReactionType.RH_POS))
        if w > 0 > lb:
            rxn_data.append(_RxnVar(i, id, lb, ub, w, _ReactionType.RH_NEG))
        if w < 0 < abs(lb) + abs(ub):
            rxn_data.append(_RxnVar(i, id, lb, ub, w, _ReactionType.RL))
    return rxn_data


class BaseModel(ABC):
    """Base class for building LP/MIP optimization models
    using a high-level API for metabolic problems.
    
    It implements the chainable methods to set-up a LP/MIP 
    metabolic network problem. Two implementations are available: 
    PicosModel and PythonMipModel. The PicosModel uses the 
    [PICOS library](https://picos-api.gitlab.io/picos/) 
    as a backend to interact with different solvers. The PythonMipModel uses 
    the [Python-MIP](https://www.python-mip.com/) library to solve the
    model using the CBC or GUROBI solvers.
    
    !!! note
        Do not try to instantiate this class directly. Use the [load()][miom.miom.load] 
        function instead. The method automatically selects the right implementation 
        depending on the solver.
        
    """
    def __init__(self, previous_step_model=None, miom_network=None, solver_name=None):
        self.previous_step_model = previous_step_model
        self.network = miom_network
        self.variables = None
        self.objective = None
        self._last_start_time = None
        self.last_solver_time = None
        if previous_step_model is not None:
            self._options = previous_step_model._options
            if miom_network is None:
                self.network = previous_step_model.network
            if solver_name is not None:
                self._options["solver"] = solver_name
        else:
            # Default recommended options
            self._options = {
                'int_tol': 1e-8,
                'feas_tol': 1e-8,
                'opt_tol': 1e-5,
                'verbosity': 0,
                'solver': solver_name
            }
        self.problem = self.initialize_problem()
        self.setup(**self._options)

    @abstractmethod
    def initialize_problem(self):
        pass

    @abstractmethod
    def get_solver_status(self):
        pass

    @property
    def status(self):
        return self.get_solver_status()
        

    @_composable
    def setup(self, **kwargs):
        """Provide the options for the solver.

        Attributes:
            int_tol (float): Integrality tolerance for integer variables.
                Defaults to 1e-8.
            feas_tol (float): Feasibility tolerance. Defaults to 1e-8.
            opt_tol (float): Relative MIP gap tolerance for MIP problems.
                Defaults to 1e-5.
            verbosity (int): Values above 0 force the backends to be verbose.
                Use a value of 1 to show useful information about the search.
                Defaults to 0.

        Returns:
            BaseModel: the configured instance of the BaseModel
        """
        if self.problem is None:
            warnings.warn("Problem cannot be configured since it was not initialized")
            return False
        options = dict()
        int_tol = kwargs["int_tol"] if "int_tol" in kwargs else None
        feas_tol = kwargs["feas_tol"] if "feas_tol" in kwargs else None
        opt_tol = kwargs["opt_tol"] if "opt_tol" in kwargs else None
        verbosity = kwargs["verbosity"] if "verbosity" in kwargs else None
        solver = kwargs["solver"] if "solver" in kwargs else None
        if int_tol is not None:
            options["int_tol"] = int_tol
        if feas_tol is not None:
            options["feas_tol"] = feas_tol
        if opt_tol is not None:
            options["opt_tol"] = opt_tol
        if verbosity is not None:
            options["verbosity"] = verbosity
        if solver is not None:
            options["solver"] = solver
        self._options.update(options)
        # Calculate a min eps value to avoid numerical issues
        self._options["_min_eps"] = np.sqrt(2*self._options["int_tol"])
        return self._options

    @_composable
    def steady_state(self):
        """Add the required constraints for finding steady-state fluxes

        The method adds the $S * V = 0$ set of constraints, where $S$
        is the stoichiometric matrix and $V$ the flux variables.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        pass

    @_composable
    def keep(self, reactions):
        """Force the inclusion of a list of reactions in the solution.

        Reactions have to be associated with positive weights in order to
        keep them in the final solution. Note that once keep() is called,
        the weights associated to the reactions selected to be kept will
        not be taken into account, as they will be forced to be kept in
        the solution.

        Args:
            reactions (list): List of reaction names, a binary vector
                indicating the reactions to keep, or a list of indexes
                with the reactions to keep.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        if self.variables.indicators is None:
            raise ValueError("""No indicator variables for reactions, 
                             transform it to a subset selection problem calling
                             subset_selection first, providing positive weights for the
                             reactions you want to keep.""")
        if reactions is None or len(reactions) == 0:
            return False
        if isinstance(reactions, str):
            reactions = [reactions]
        # Check if it's a binary vector
        if len(reactions) == self.network.num_reactions and max(reactions)==1:
            reactions = np.where(reactions==1)[0]
        # Check if the reactions have an indicator variable
        reactions = set(self.network.find_reaction(rxn)[0] for rxn in reactions)
        available = set(rxn.index for rxn in self.variables.assigned_reactions if rxn.cost > 0)
        diff = reactions - available
        if len(diff) != 0:
            raise ValueError(f"""Only reactions associated with positive weights
                             can be forced to be selected. The following reaction 
                             indexes have no indicator variables or are associated with 
                             negative weights: {diff}.""")
        valid_rxn_idx = reactions & available
        # Get the indexes of the indicator variables
        idxs = [i for i, r in enumerate(self.variables.assigned_reactions)
                if r.index in valid_rxn_idx]
        return dict(idxs=idxs, valid_rxn_idx=valid_rxn_idx)

    @_composable
    def subset_selection(self, rxn_weights, direction='max', eps=1e-2, strengthen=True):
        """Transform the current model into a subset selection problem.

        The subset selection problem consists of selecting the positive weighted
        reactions and remove the negative weighted reactions, subject to steady
        state constraints and/or additional constraints on fluxes, and maximizing
        the weighted sum of the (absolute) weights for the successfully selected reactions
        (with positive weights) and the successfully removed reactions (with negative
        weights). Selected reactions are forced to have an absolute flux value greater 
        or equal to the threshold `eps` (1e-2 by default). Removed reactions should have a
        flux equal to 0.

        Each reaction is associated with a weight (positive or negative) provided
        in the parameter `rxn_weights`, and the objective is to select the reactions 
        that optimizes the following expression:

        $$
        f(x) = \sum_i^n |w_i| * x_i
        $$

        where $x_i$ are the indicator variables for the reactions $i$ and $w_i$ are
        the weights for the reactions associated to the indicator variable. Indicator
        variables are automatically created for each reaction associated to a non-zero
        weight. Two (mutually exclusive) indicator variables are used for positive weighted 
        reactions that are reversible to indicate whether there is positive or negative flux 
        through the reaction. A single indicator variable is created for positive weighted
        non-reversible reactions, to indicate if the reaction is selected (has a non-zero 
        flux greater or equal to `eps`) in which case the indicator variable is 1, 
        or 0 otherwise. 
        
        A single binary indicator variable is also created for negative weighted reactions, 
        indicating whether the reaction was not selected (i.e, has 0 flux, in which case the
        indicator variable is 1) or not (in which case the indicator variable is 0).

        Args:
            rxn_weights (list): List of weights for each reaction. If a single value
                is provided, it is assumed to be the weight for all reactions.
            eps (float, optional): Min absolute flux value for weighted reactions
                to consider them active or inactive. Defaults to 1e-2.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        # Calculate min valid EPS based on integrality tolerance
        min_eps = self._options["_min_eps"]
        if eps < min_eps:
            warnings.warn(f"Minimum epsilon value below min. allowed value, changed to {min_eps}.")
        eps = max(eps, min_eps)
        if not isinstance(rxn_weights, Iterable):
            rxn_weights = [rxn_weights] * self.network.num_reactions
        rxnw = _weighted_rxns(self.network.R, rxn_weights)
        if self.variables.indicators is None:
            self.variables._assigned_reactions = rxnw
            return dict(eps=eps, direction=direction, strengthen=strengthen)
        else:
            raise ValueError("The current model is already a subset selection problem.")

    @_composable
    def set_flux_bounds(self, rxn_id, min_flux=None, max_flux=None):
        """Change the flux bounds of a reaction.

        Args:
            rxn_id (str/int): name or id of the reaction to change
            min_flux (float, optional): Min flux value. Defaults to None.
            max_flux (float, optional): Max flux value. Defaults to None.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        i, _ = self.network.find_reaction(rxn_id)
        return i

    @_composable
    def add_constraints(self, constraints):
        """Add a list of constraint to the model

        Args:
            constraints (list): List of expressions with the
                constraints.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        for c in constraints:
            self.add_constraint(c)
        return len(constraints) > 0

    @_composable
    def exclude(self, indicator_values=None):
        """Exclude a solution from the set of feasible solutions.

        If the problem is a subset selection problem, it adds a new constraint
        to exclude the given values (or the current values of the indicator variables)
        from the set of feasible solutions. This is useful to force the solver to find
        alternative solutions in a manual fashion. https://doi.org/10.1371/journal.pcbi.1008730
        

        Args:
            indicator_values (list, optional): List of values for each indicator variable. Defaults to None.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        if self.variables.indicators is None:
            raise ValueError("""The optimization model does not contain indicator variables.
            Make sure that the problem is a subset selection problem.""")
        if indicator_values is None:
            indicator_values = np.array(self.variables.indicator_values)
        return dict(indicator_values=indicator_values)

    @_composable
    def add_constraint(self, constraint):
        """Add a specific constraint to the model.
        The constraint should use existing variables already included in the model.

        Args:
            constraint: affine expression using model's variables.
        """
        pass

    @_composable
    def set_objective(self, cost_vector, variables, direction='max'):
        """Set the optmization objective of the model.

        Args:
            cost_vector (Iterable): List with the cost for each variable
            variables (Iterable): Variables used for the objective function
            direction (str, optional): Optimization direction (min or max). Defaults to 'max'.
        """
        if self.objective is not None:
            warnings.warn("Previous objective changed")
        self.objective = (cost_vector, variables, direction)

    def set_rxn_objective(self, rxn, direction='max'):
        """Set a flux objective

        Maximize or minimize the flux through a given reaction.

        Args:
            rxn (str): Name of the reaction to use for the optimization
            direction (str, optional): Minimization or maximization of 
                the flux ('min' or 'max'). Defaults to 'max'.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        i = self.network.get_reaction_id(rxn)
        cost = np.zeros(self.network.num_reactions)
        cost[i] = 1
        self.set_objective(cost, self.variables.fluxvars, direction=direction)
        return self


    def set_fluxes_for(self, reactions, tolerance=1e-6):
        warnings.warn("This method was renamed to fix_fluxes. It will dissappear in the v0.9.0", DeprecationWarning)
        return self.fix_fluxes(reactions, tolerance)


    def fix_fluxes(self, reactions, tolerance=1e-6):
        """Force the flux of certain reactions to match current values.

        After calling `.solve()` for a flux optimization problem (e.g. FBA), this
        method adds a new constraint to force the flux of the given reactions to
        match the current flux values found after the optimization.

        This is interesting for example to implement methods like sparse-FBA, where
        the optimization is no longer the flux but the number of active reactions,
        and a new constraint is required to preserve optimality of fluxes.

        Args:
            reactions (list): reaction or list of reactions
            tolerance (float, optional): Tolerance for the flux values
                (a solution is valid if the flux is within optimal value +/- tol. 
                Defaults to 1e-6.

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        if isinstance(reactions, str):
            reactions = [reactions]

        for rid in reactions:
            if isinstance(rid, np.ndarray):
                # TODO: Test this path
                rid = rid['id']
            idx, rxn = self.network.find_reaction(rid)
            lb = max(rxn['lb'], self.variables.flux_values[idx] - tolerance)
            ub = min(rxn['ub'], self.variables.flux_values[idx] + tolerance)
            self.add_constraint(self.variables.fluxvars[idx] >= lb)
            self.add_constraint(self.variables.fluxvars[idx] <= ub)
        return self

    @_composable
    def reset(self):
        """Resets the original problem (removes all modifications)

        Returns:
            BaseModel: instance of BaseModel with the modifications applied.
        """
        if self.problem is None:
            warnings.warn("Problem is not initialized, nothing to reset")
            return False
        else:
            return True

    def obtain_subnetwork(
            self,
            extraction_mode=ExtractionMode.ABSOLUTE_FLUX_VALUE,
            comparator=Comparator.GREATER_OR_EQUAL,
            value=1e-8
    ):
        """Same as [select_subnetwork][miom.miom.BaseModel] but returns the network instead.

        See [select_subnetwork][miom.miom.BaseModel] for a detailed description of the method.

        Returns:
            MiomNetwork: A MiomNetwork with the selected subnetwork.
        """
        # If indicators are present and assigned,
        # take the subset of the network for which
        # the indicators of positive weighted reactions
        # are equal to 1
        if extraction_mode == ExtractionMode.ABSOLUTE_FLUX_VALUE:
            variables = self.variables.flux_values
            if variables is None:
                raise ValueError("The model does not contain flux variables. "
                                 "You need to call first steady_state() to add "
                                 "the flux variables")
        elif extraction_mode == ExtractionMode.INDICATOR_VALUE:
            variables = self.variables.indicator_values
            if variables is None:
                raise ValueError("The model does not contain indicator variables. "
                                 "You need to transform it to a subset selection problem "
                                 "by invoking subset_selection() first.")
        elif extraction_mode == ExtractionMode.REACTION_ACTIVITY:
            variables = self.variables.reaction_activity
            if variables is None:
                raise ValueError("The model does not contain reaction activity variables. "
                                 "You need to transform it to a subset selection problem "
                                 "by invoking subset_selection() first.")
        else:
            raise ValueError("Invalid extraction mode")

        if comparator == Comparator.GREATER_OR_EQUAL:
            selected = np.where(np.abs(variables) >= value)[0]
        elif comparator == Comparator.LESS_OR_EQUAL:
            selected = np.where(np.abs(variables) <= value)[0]
        else:
            raise ValueError("Invalid comparison")
        # TODO: Assigned reactions work only for indicator variables
        if extraction_mode == ExtractionMode.INDICATOR_VALUE:
            rxns = [self.variables.assigned_reactions[x] for x in selected]
            selected_idx = [rxn.index for rxn in rxns]
        elif extraction_mode == ExtractionMode.ABSOLUTE_FLUX_VALUE:
            selected_idx = selected
        else:
            raise NotImplementedError("Only indicator variables and absolute flux values are supported")
        return self.network.subnet(selected_idx)

    @_composable
    def select_subnetwork(
            self,
            mode=ExtractionMode.ABSOLUTE_FLUX_VALUE,
            comparator=Comparator.GREATER_OR_EQUAL,
            value=1e-8
    ):
        """Select a subnetwork and create a new BaseModel to operate on it.

        The new instance of the BaseModel is a new problem instance with no constraints.
        If the idea is to perform FBA simulations on this new subnetwork, remember to add
        the new constraints, especially the `steady_state`.

        Args:
            mode (ExtractionMode, optional): Method used to extract the subnetwork 
                (based on flux values or using the indicator values). 
                Defaults to ExtractionMode.ABSOLUTE_FLUX_VALUE.
            comparator (Comparator, optional): Comparator for the selected mode. 
                Defaults to Comparator.GREATER_OR_EQUAL.
            value (float, optional): Value threshold for the mode and comparator selected. 
                Defaults to 1e-8.

        Returns:
            [type]: [description]
        """
        return self.obtain_subnetwork(extraction_mode=mode,
                                      comparator=comparator,
                                      value=value)

    def get_values(self):
        """Get the values for the variables

        Returns:
            tuple: (V, X) where V are the flux values and X are the indicator values
                (if the problem is a MIP problem, for example if `subset_selection` was
                called)
        """
        return self.variables.values()
        

    def get_fluxes(self, reactions=None):
        """Get the flux values.

        Args:
            reactions (list, optional): Reaction or subset of reactions
                For which to obtain the flux values. Defaults to None.

        Returns:
            list: List with the flux values for all or the selected reactions.
        """
        if isinstance(reactions, str):
            return self.variables.flux_values[self.network.get_reaction_id(reactions)]
        if isinstance(reactions, Iterable):
            return {
                r['id']: self.variables.flux_values[self.network.get_reaction_id(r['id'])]
                for r in reactions
            }
        if reactions is None:
            return self.variables.flux_values
        else:
            raise ValueError("reactions should be an iterable of strings or a single string")

    @_composable
    def solve(self, verbosity=None, max_seconds=None):
        """Solve the current model and assign the values to the variables of the model.

        Args:
            verbosity (int, optional): Level of verbosity for the solver. 
                Values above 0 will force the backend to show output information of the search. Defaults to None.
            max_seconds (int, optional): Max time in seconds for the search. Defaults to None.
        """
        self._last_start_time = perf_counter()

    @_composable
    def copy(self):
        pass


class PicosModel(BaseModel):
    def __init__(self, previous_step_model=None, miom_network=None, solver_name=None):
        super().__init__(previous_step_model=previous_step_model,
                         miom_network=miom_network,
                         solver_name=solver_name)
        self.variables = _PicosVariables()

    def initialize_problem(self):
        return pc.Problem()

    def _setup(self, *args, **kwargs):
        self.problem.options["verbosity"] = kwargs["verbosity"] if "verbosity" in kwargs else self._options["verbosity"]
        self.problem.options["rel_bnb_opt_tol"] = kwargs["opt_tol"] if "opt_tol" in kwargs else self._options["opt_tol"]
        self.problem.options["integrality_tol"] = kwargs["int_tol"] if "int_tol" in kwargs else self._options["int_tol"]
        self.problem.options["abs_prim_fsb_tol"] = kwargs["feas_tol"] if "feas_tol" in kwargs else self._options["feas_tol"]
        self.problem.options["abs_dual_fsb_tol"] = self.problem.options["abs_prim_fsb_tol"]
        self.problem.options["solver"] = kwargs["solver"] if "solver" in kwargs else self._options["solver"]
        return True

    def _reset(self, **kwargs):
        self.problem.reset()
        return True

    def _keep(self, *args, **kwargs):
        idx = kwargs["idxs"]
        valid_rxn_idx = kwargs["valid_rxn_idx"]
        n = len(valid_rxn_idx)
        C = pc.Constant("C", value=[1] * len(idx))
        # The sum should be equal to the number of different reactions
        self.add_constraint(C.T * self.variables.indicators[idx] >= n)

    def _exclude(self, *args, **kwargs):
        values = np.round(kwargs["indicator_values"])
        b = sum(values) - 1
        idx_one = np.flatnonzero(values==1)
        idx_zero = np.flatnonzero(values==0)
        s = sum(self.variables.indicators[idx_one] // self.variables.indicators[idx_zero])
        self.add_constraint(s <= b)

    def _set_flux_bounds(self, *args, **kwargs):
        i = kwargs["_parent_result"]
        min_flux = kwargs["min_flux"] if "min_flux" in kwargs else None
        max_flux = kwargs["max_flux"] if "max_flux" in kwargs else None
        if min_flux is not None:
            self.variables.fluxvars._lower[i] = min_flux
        if max_flux is not None:
            self.variables.fluxvars._upper[i] = max_flux

    def _subset_selection(self, *args, **kwargs):
        # Convert to a weighted subset selection problem
        eps = kwargs["eps"]
        weighted_rxns = self.variables.assigned_reactions
        P = self.problem
        # Build the MIP problem
        V = self.variables.fluxvars
        X = pc.BinaryVariable("X", shape=(len(weighted_rxns),))
        C = pc.Constant("C", value=np.array([abs(rxn.cost) for rxn in weighted_rxns]))
        # NOTE: For some reason, PICOS have issues parsing the priority operators, so constraints
        # are split into many parts, otherwise it does not work. The matrix notation is more
        # efficient and takes less time and memory than using the add_list_of_constraints method in a loop.
        var_data = _get_rxn_var_data(weighted_rxns, _ReactionType.RH_POS)
        if var_data is not None:
            I, J, LB, UB = var_data
            EA = eps - LB
            EB = X[I]^EA
            expr = V[J] >= EB + LB
            P.add_constraint(expr)
        var_data = _get_rxn_var_data(weighted_rxns, _ReactionType.RH_NEG)
        if var_data is not None:
            I, J, LB, UB = var_data
            EA = eps + UB
            EB = X[I]^EA
            expr = V[J] <= UB - EB
            P.add_constraint(expr)
        var_data = _get_rxn_var_data(weighted_rxns, _ReactionType.RL)
        if var_data is not None:
            I, J, LB, UB = var_data
            EA = 1 - X[I]
            EB = EA^UB
            EC = EA^LB
            P.add_constraint(EB - V[J] >= 0)
            P.add_constraint(EC - V[J] <= 0)
        # Strengthen the formulation
        if kwargs["strengthen"] is True:
            idx = [np.array(x) for x in list(zip(*_get_reversible_vars(weighted_rxns)))]
            if len(idx) > 0:
                I, J = idx
                P.add_constraint(X[I] + X[J] <= 1)
        P.set_objective(kwargs['direction'], C.T * X)
        self.variables._indicator_vars = X
        return True

    def _solve(self, **kwargs):
        max_seconds = kwargs["max_seconds"] if "max_seconds" in kwargs else None
        verbosity = kwargs["verbosity"] if "verbosity" in kwargs else None
        init_max_seconds = self.problem.options["timelimit"]
        init_verbosity = self.problem.options["verbosity"]
        if max_seconds is not None:
            self.problem.options["timelimit"] = max_seconds
        if verbosity is not None:
            self.problem.options["verbosity"] = verbosity
        self.solutions = self.problem.solve()
        self.problem.options["timelimit"] = init_max_seconds
        self.problem.options["verbosity"] = init_verbosity
        self.last_solver_time = perf_counter() - self._last_start_time
        return True

    def _add_constraint(self, constraint, **kwargs):
        self.problem.add_constraint(constraint)

    def _steady_state(self, **kwargs):
        S = pc.Constant("S", self.network.S)
        V = pc.RealVariable("V", (self.network.S.shape[1], 1),
                            lower=[rxn['lb'] for rxn in self.network.R],
                            upper=[rxn['ub'] for rxn in self.network.R])
        self.add_constraint(S * V == 0)
        self.variables._flux_vars = V
        return True

    def _set_objective(self, cost_vector, variables, **kwargs):
        if "direction" in kwargs:
            direction = kwargs["direction"]
        else:
            direction = "max"
        C = pc.Constant("C", value=cost_vector)
        self.problem.set_objective(direction, C.T * variables)
        return True

    def _select_subnetwork(self, **kwargs):
        m = kwargs["_parent_result"]
        return PicosModel(previous_step_model=self, miom_network=m)

    def _copy(self, **kwargs):
        m = PicosModel(previous_step_model=self)
        m.problem = self.problem.copy()
        variables = m.problem._mutables.keys()
        m.variables._assigned_reactions = self.variables._assigned_reactions
        if "X" in variables:
            m.variables._indicator_vars = m.problem._mutables.get("X")
        if "V" in variables:
            m.variables._flux_vars = m.problem._mutables.get("V")
        return m

    def get_solver_status(self):
        # Check if attribute solutions exist
        if not hasattr(self, "solutions"):
            status = "empty"
        else:
            status = self.solutions.claimedStatus
        return {
            "status": status,
            "objective_value": self.problem.value,
            "elapsed_seconds": self.last_solver_time
        }


class PythonMipModel(BaseModel):
    
    _status_eq = {
        mip.OptimizationStatus.OPTIMAL: Status.OPTIMAL,
        mip.OptimizationStatus.FEASIBLE: Status.FEASIBLE,
        mip.OptimizationStatus.INFEASIBLE: Status.INFEASIBLE,
        mip.OptimizationStatus.UNBOUNDED: Status.UNBOUNDED,
        mip.OptimizationStatus.INT_INFEASIBLE: Status.INFEASIBLE,
        mip.OptimizationStatus.NO_SOLUTION_FOUND: Status.ILLPOSED,
        mip.OptimizationStatus.LOADED: Status.EMPTY,
        mip.OptimizationStatus.CUTOFF: Status.PREMATURE
    }

    def __init__(self, previous_step_model=None, miom_network=None, solver_name=None):
        super().__init__(previous_step_model=previous_step_model,
                         miom_network=miom_network,
                         solver_name=solver_name)
        self.variables = _PythonMipVariables()

    def initialize_problem(self):
        solver = self._options["solver"]
        if solver is None:
            solver = 'cbc'
        else:
            solver = solver.lower()
        if solver == 'gurobi':
            mip_solver = mip.GRB
        elif solver == 'cbc':
            mip_solver = mip.CBC
        else:
            raise ValueError("Only gurobi and cbc are supported with Python-MIP")
        return mip.Model("model", solver_name=mip_solver)

    def _setup(self, *args, **kwargs):
        self.problem.max_mip_gap = kwargs["opt_tol"] if "opt_tol" in kwargs else self._options["opt_tol"]
        self.problem.max_gap = self.problem.max_mip_gap
        self.problem.infeas_tol = kwargs["feas_tol"] if "feas_tol" in kwargs else self._options["feas_tol"]
        self.problem.opt_tol = self.problem.infeas_tol
        self.problem.integer_tol = kwargs["int_tol"] if "int_tol" in kwargs else self._options["int_tol"]
        self.problem.verbose = kwargs["verbosity"] if "verbosity" in kwargs else self._options["verbosity"]
        self.problem.store_search_progress_log = True
        self.problem.threads = -1
        return True

    def _reset(self, **kwargs):
        self.proble.clear()
        return True

    def _keep(self, *args, **kwargs):
        idx = kwargs["idxs"]
        valid_rxn_idx = kwargs["valid_rxn_idx"]
        n = len(valid_rxn_idx)
        variables = [self.variables.indicators[i] for i in idx]
        # The sum should be equal to the number of different reactions
        # Note that if a reaction is blocked and in the core, the 
        # resulting problem is infeasible
        self.add_constraint(
            mip.xsum((v for v in variables)) >= n
        )
        return True 

    def _exclude(self, *args, **kwargs):
        values = kwargs["indicator_values"]
        b = sum(values) - 1
        x1 = mip.xsum(v for i, v in enumerate(self.variables.indicators) if values[i] == 1)
        x2 = mip.xsum(v for i, v in enumerate(self.variables.indicators) if values[i] == 0)
        self.add_constraint(x1 + x2 <= b)
        
    def _set_flux_bounds(self, *args, **kwargs):
        i = kwargs["_parent_result"]
        min_flux = kwargs["min_flux"] if "min_flux" in kwargs else None
        max_flux = kwargs["max_flux"] if "max_flux" in kwargs else None
        if min_flux is not None:
            self.variables.fluxvars[i].lb = min_flux
        if max_flux is not None:
            self.variables.fluxvars[i].ub = max_flux

    def _subset_selection(self, *args, **kwargs):
        eps = kwargs["eps"]
        direction = kwargs["direction"]
        weighted_rxns = self.variables.assigned_reactions
        P = self.problem
        V = self.variables.fluxvars
        X = [self.problem.add_var(var_type=mip.BINARY, 
            name=f"X_{i}_rxn_{rxn.index}") for i, rxn in enumerate(weighted_rxns)]
        C = [abs(rxn.cost) for rxn in weighted_rxns]
        for i, rd in enumerate(weighted_rxns):
            if rd.type is _ReactionType.RH_POS:
                self.add_constraint(V[rd.index] >= X[i] * (eps - rd.lb) + rd.lb)
            elif rd.type is _ReactionType.RH_NEG:
                self.add_constraint(V[rd.index] <= rd.ub - X[i] * (rd.ub + eps))
            elif rd.type is _ReactionType.RL:
                self.add_constraint((1 - X[i]) * rd.ub - V[rd.index] >= 0)
                self.add_constraint((1 - X[i]) * rd.lb - V[rd.index] <= 0)
            else:
                raise ValueError("Unknown type of reaction")
        # Strengthen the problem for RH+ and RH-
        if kwargs["strengthen"] is True:
            for i, j in _get_reversible_vars(weighted_rxns):
                self.add_constraint(X[i] + X[j] <= 1)
        P.objective = mip.xsum((c * X[i] for i, c in enumerate(C)))
        if direction == "max" or direction == "maximize":
            P.sense = mip.MAXIMIZE
        elif direction == "min" or direction == "minimize":
            P.sense = mip.MINIMIZE
        else:
            raise ValueError("Unknown optimization sense")
        self.variables._indicator_vars = X
        return True

    def _solve(self, **kwargs):
        max_seconds = kwargs["max_seconds"] if "max_seconds" in kwargs else 10 * 60
        verbosity = kwargs["verbosity"] if "verbosity" in kwargs else None
        init_verbosity = self.problem.verbose
        if verbosity is not None:
            self.setup(verbosity=verbosity)
        solutions = self.problem.optimize(max_seconds=max_seconds)
        self.setup(verbosity=init_verbosity)
        self.last_solver_time = perf_counter() - self._last_start_time
        return True

    def _add_constraint(self, constraint, **kwargs):
        self.problem += constraint

    def _steady_state(self, **kwargs):
        V = [self.problem.add_var(lb=rxn['lb'], ub=rxn['ub'], name=f"V_{i}") for i, rxn in enumerate(self.network.R)]
        for i in range(self.network.S.shape[0]):
            non_zero_cols = np.flatnonzero(self.network.S[i, :] != 0)
            if len(non_zero_cols) > 0:
                self.problem += mip.xsum(self.network.S[i, j] * V[j]
                    for j in non_zero_cols) == 0
        self.variables._flux_vars = V
        return True

    def _set_objective(self, cost_vector, variables, **kwargs):
        if "direction" in kwargs:
            direction = kwargs["direction"]
        else:
            direction = "max"
        if direction == "max":
            self.problem.sense = mip.MAXIMIZE
        else:
            self.problem.sense = mip.MINIMIZE
        self.problem.objective = (
            mip.xsum(
                float(cost_vector[i]) * variables[i] for i in range(len(variables))
            ) 
        )      
        return True

    def _select_subnetwork(self, **kwargs):
        return PythonMipModel(previous_step_model=self, miom_network=kwargs["_parent_result"])

    def _copy(self, **kwargs):
        m = PythonMipModel(previous_step_model=self)
        m.problem = self.problem.copy()
        m.variables._assigned_reactions = self.variables._assigned_reactions
        fluxvars = [v for v in m.problem.vars if v.name.startswith("V_")]
        indicators = [v for v in m.problem.vars if v.name.startswith("X_")]
        if len(fluxvars) > 0:
            m.variables._flux_vars = fluxvars
        if len(indicators) > 0:
            m.variables._indicator_vars = indicators
        return m
        

    def get_solver_status(self):
        return {
            "status": PythonMipModel._status_eq[self.problem.status] \
                if self.problem.status in PythonMipModel._status_eq \
                     else Status.UNKNOWN,
            "objective_value": self.problem.objective_value,
            "elapsed_seconds": self.last_solver_time
        }

