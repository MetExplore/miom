from miom import miom, load_gem, Solvers

# Load a genome-scale metabolic network. You can load SMBL or Matlab metabolic networks
# as well using the same method, but it requires to have the cobratoolbox python library
# installed.
network = load_gem("https://github.com/pablormier/miom-gems/raw/main/gems/mus_musculus_iMM1865.miom")
target_rxn = "BIOMASS_reaction"
# Create the optimization problem with miom and solve
model = (miom(network)
        .steady_state()
        .set_rxn_objective(target_rxn)
        .solve(verbosity=1))
print("Optimal flux:", model.get_fluxes(target_rxn), "mmol/(h·gDW)")
# Show reactions with non-zero flux
V, _ = model.get_values()
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

print("Number of reactions with non-zero flux:", sum(abs(V) > 1e-8))
