# Sparse FBA (MIP)

The sparse FBA problem consists of finding the optimal flux value for a reaction minimizing the number of reactions with non-zero flux (l0-norm sparsity). Note that minimizing the l0-norm is
a NP-hard problem, and so obtaing an optimal solution is not possible in many cases. The COBRA
Toolbox includes different LP heuristics to minimize different approximations of the l0-norm. However, using a MIO formulation, it's is possible to obtain a solution which is close to optimal, with a lower number of reactions than the LP approximations, and very quickly.

Here is the implementation of sparse-FBA with MIOM:

```python
import miom

# Load a genome-scale metabolic network. You can load SMBL or Matlab metabolic networks
# as well using the same method, but it requires to have the cobratoolbox python library
# installed.
network = miom.mio.load_gem("https://github.com/pablormier/miom-gems/raw/main/gems/mus_musculus_iMM1865.miom")

# Create the sparse FBA problem to get a solution that maximizes
# the optimal flux through the BIOMASS_reaction minimizing the
# number of active reactions. The solution should be not more than
# 5% of the optimal solution (opt_tol = 0.05).
V, X = (miom
        .load(network, solver=miom.Solvers.GUROBI_PYMIP)
        # Set-up the solver options
        .setup(int_tol=1e-8, opt_tol=0.05, verbosity=1)
        # Add the steady-state constraints (S*V = 0)
        .steady_state()
        # Set the reaction to optimize using FBA
        .set_rxn_objective('BIOMASS_reaction')
        # Solve the FBA (LP) problem (optimal flux)
        .solve()
        # Add a constraint to force a flux through
        # the reaction equal to the optimal flux
        .set_fluxes_for('BIOMASS_reaction')
        # Convert to a MIO problem (best subset selection)
        # Each reaction in the network is associated
        # with a negative weight of -1. The optimization
        # problem now tries to minimize the selection of
        # reactions with negative weights (respecting
        # the previous constraints). Since each reaction 
        # has a weight of -1, all reactions are equally 
        # important in the optimization problem.
        # This is exactly the optimal sparse-FBA problem
        # with l0-norm: minimize the number of reactions
        # but mantaining the optimal FBA flux possible.
        .subset_selection(-1)
        # Solve the MIO problem (note that it's NP-Hard)
        # You can control the optimality by changing the
        # opt_tol in the .setup() method.
        .solve()
        # Return the flux values (V) and the binary 
        # indicator values (X)
        .get_values())

# Show reactions with a flux > 1e-7
print("Number of reactions with flux above +/- 1e-8:", sum(abs(V) > 1e-8))

# Count reactions with an indicator value of 0 (active). Note that since
# the weights of the reactions are negative (for all rxns), an indicator
# value of 1 corresponds to a succesful removed reaction not included in
# the solution, and a value of 0 to a reaction that could not be removed.
print("Indicator variables with a value of 1 (selected rxns):", sum(X < 0.5))

# Show the flux value of the biomass reaction
print("Optimal biomass flux:", V[m.get_reaction_id("BIOMASS_reaction")], "mmol/(h·gDW)")

```
