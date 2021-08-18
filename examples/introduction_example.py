import miom

# Load a genome-scale metabolic network. You can load SMBL or Matlab metabolic networks
# as well using the same method, but it requires to have the cobratoolbox python library
# installed.
network = miom.load_gem("@mus_musculus_iMM1865.miom")
target_rxn = "BIOMASS_reaction"
# Create the optimization problem with miom and solve
model = (miom
        .load(network)
        .steady_state()
        .set_rxn_objective(target_rxn)
        .solve(verbosity=1))
print("Optimal flux:", model.get_fluxes(target_rxn), "mmol/(h·gDW)")

# Show reactions with non-zero flux
V = model.get_fluxes()
print("Number of reactions with flux above +/- 1e-8:", sum(abs(V) > 1e-8))

V, X = (model
        # Set the MIP Gap tolerance to 5%, using the default solver
        # (COIN-OR CBC included with the Python-MIP lib)
        .setup(opt_tol=0.05, verbosity=1)
        # Set the fluxes to the optimal value found
        .set_fluxes_for('BIOMASS_reaction')
        # Transform to a subset selection problem
        # adding a negative weight to all reactions
        # (to remove them from the final solution)
        .subset_selection(-1)
        # Solve the MIO problem
        .solve()
        # Get continuos vars (fluxes) and binary vars
        .get_values())

print("Number of reactions with an absolute flux value above 1e-8:", sum(abs(V) > 1e-8))
print("Active reactions:", sum(1 if activity == 1 else 0 for activity in model.variables.reaction_activity))
print("Inconsistencies:", sum(1 if activity != activity else 0 for activity in model.variables.reaction_activity))
print("Solver status:", model.get_solver_status())