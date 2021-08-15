# Exact Fastcore (MIP)

Example of implementation of the Fastcore algorithm with MIOM to extract a context specific model.
Fastcore defines a set of core reactions that is forced to be active (carry a non-zero flux in steady state conditions), and minimizes the set of non-core reactions.

!!! cite
    Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact context-specific metbolic network models. [PLoS Comput. Biol. 10, e1003424](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424).


The Fastcore algorithm is greedy approximation of the exact problem which can be modelled as a MIP problem. With MIOM, the exact problem can be defined in a few lines, using the `opt_tol` parameter to control the level of optimality required:

```python
import miom
import numpy as np

# Use the flux-consistent subnetwork (fcm) of the Human1 GEM model 
m = miom.mio.load_gem('https://github.com/pablormier/miom-gems/raw/main/gems/homo_sapiens_human1_fcm.miom')
# Select reactions from the cholesterol metabolism as the core reactions to keep
core_rxn = m.find_reactions_from_pathway("Cholesterol metabolism")
print(sum(core_rxn))

# Assign a negative weight for reactions not in the core
weights = -1 * np.ones(m.num_reactions)
weights[core_rxn == 1] = 1

# Exact-Fastcore
fmc = (miom
        .load(m, solver=miom.Solvers.GUROBI_PYMIP)
        .setup(opt_tol=0.01)
        .steady_state()
        .subset_selection(weights)
        .keep(core_rxn == 1)
        .solve(verbosity=1)
        .select_subnetwork(
            mode=miom.ExtractionMode.ABSOLUTE_FLUX_VALUE,
            comparator=miom.Comparator.GREATER_OR_EQUAL,
            value=1e-8
        )
        .network
)
print(fmc.num_reactions)
```

Note that the algorithm requires that the metabolic network used is flux consistent. If blocked reactions are included in the core set, the MIP becomes infeasible, as those reactions cannot be selected with a non-zero flux in steady state conditions.

MIOM includes an implementation of the [swiftcore algorithm](https://mtefagh.github.io/swiftcore/) to obtain flux consistent metabolic networks. 