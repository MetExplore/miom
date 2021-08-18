import miom
import numpy as np

# Load a model
m = miom.load_gem('@homo_sapiens_human1_fcm.miom')
print("Reactions in the consistent iHuman GEM:", m.num_reactions)
core_rxn = m.find_reactions_from_pathway("Cholesterol metabolism")
print("Num. of core reactions:", sum(core_rxn))

# Assign a negative weight for reactions not in the core and a positive value
# for reactions in the core.
weights = -1 * np.ones(m.num_reactions)
weights[core_rxn == 1] = 1

# Weighted Exact Fastcore algorithm with MIOM:
fmc = (miom.load(m, solver=miom.Solvers.COIN_OR_CBC)
          .setup(opt_tol=0.05)
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
print("Number of reactions in the fastcore network:", fmc.num_reactions)
print("Number of reactions selected from the core:", sum(fmc.find_reactions_from_pathway("Cholesterol metabolism")))