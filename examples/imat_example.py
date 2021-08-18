import miom

# Example implementation of iMAT using MIOM.
# Note that this implementation supports also custom weights for the reactions
# By default, weights for reactions in the highly expressed set are set to 1,
# and reactions in the lowly expressed set are set to -1, so the objective
# function is exactly the same as the original iMAT.

# Use the iHuman-GEM model
m = miom.load_gem('@homo_sapiens_human1.miom')

# Add all the reactions from the Cholesterol pathway to the highly expressed set
RH = m.find_reactions_from_pathway("Cholesterol metabolism")
# Add reactions from the pyruvate metabolism to the lowly expressed set
RL = -1 * m.find_reactions_from_pathway("Pyruvate metabolism")
w = RH + RL
print("RH:", sum(RH), "RL:", sum(abs(RL)))

m = (miom(m, solver=miom.Solvers.COIN_OR_CBC)
     .setup(int_tol=1e-8, opt_tol=0.01, verbosity=1)
     .steady_state()
     .subset_selection(w)
     .solve(max_seconds=30)
     .select_subnetwork(
          mode=miom.ExtractionMode.ABSOLUTE_FLUX_VALUE,
          comparator=miom.Comparator.GREATER_OR_EQUAL,
          value=1e-8
     )
     .network)
print("Number of reactions in the subnetwork with active flux:", m.num_reactions)

