import miom

# Load a model
m = miom.mio.load_gem('@mus_musculus_iMM1865.miom')
print("Num. of reactions in the network:", m.num_reactions)

# Solve with Gurobi, CPLEX or CBC (other MIP solvers struggle with mediudm/large networks)
V, X = (miom
        .load(network=m, solver=miom.Solvers.COIN_OR_CBC)
        # Config the solver (e.g. set the optimality tolerance for MIP problems)
        .setup(int_tol=1e-8, opt_tol=0.05)
        # Add steady-state constraints to the model (S * V = 0)
        .steady_state()
        # Set the BIOMASS_reaction as the objective flux value to maximize
        .set_rxn_objective('BIOMASS_reaction')
        # Solve the LP problem (showing info about the search process).
        # Note that verbosity can be provided in the setup() method as well.
        .solve(verbosity=1)
        # Set the flux of the biomass reaction to the optimal solution found
        .set_fluxes_for('BIOMASS_reaction')
        # Transform to a subset selection problem (MIP).
        # Assign a weight of -1 for all reactions (can be any arbitrary negative number
        # if the importance of each reaction is the same). Use different values to weight
        # more some reactions over others. The eps value indicates the min flux value 
        # to consider a reaction as active. Note that very small values can lead to nuemrical
        # issues. The min eps value depends in the int_tol parameter of the solver. A min eps
        # value is automatically calculated based on the solver's settings.
        .subset_selection(-1)
        # Solve the MIP problem
        .solve(verbosity=1)
        .get_values())

# Show reactions with a flux > 1e-8
print("Number of reactions with flux above +/- 1e-8:", sum(abs(V)>1e-8))

# Count reactions with an indicator value of 0 (active). Note that since
# the weights of the reactions are negative (for all rxns), an indicator
# value of 1 corresponds to a succesful removed reaction not included in
# the solution, and a value of 0 to a reaction that could not be removed.
print("Indicator variables with a value of 1 (selected rxns):", sum(X < 0.5))

# Show the flux value of the biomass reaction
print("Optimal biomass flux:", V[m.get_reaction_id("BIOMASS_reaction")], "mmol/(h·gDW)")
