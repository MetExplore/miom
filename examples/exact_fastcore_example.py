from miom import miom, Solvers
from miom.mio import load_gem
from miom.miom import Comparator, ExtractionMode
import numpy as np

# Load a model
m = load_gem('https://github.com/pablormier/miom-gems/raw/main/gems/homo_sapiens_human1_fcm.miom')
print("Reactions in the consistent iHuman GEM:", m.num_reactions)
core_rxn = m.find_reactions_from_pathway("Cholesterol metabolism")
print("Num. of core reactions:", sum(core_rxn))
# Assign a negative weight for reactions not in the core
weights = -1 * np.ones(m.num_reactions)
weights[core_rxn == 1] = 1

# Fastcore
fmc = (miom(m, solver=Solvers.GUROBI_PYMIP)
          .setup(opt_tol=0.05)
          .steady_state()
          .subset_selection(weights)
          .keep(np.where(core_rxn == 1)[0])
          .solve(verbosity=1)
          .select_subnetwork(
               mode=ExtractionMode.ABSOLUTE_FLUX_VALUE,
               comparator=Comparator.GREATER_OR_EQUAL,
               value=1e-8
          )
          .network
)
print("Number of reactions in the fastcore network:", fmc.num_reactions)
print("Number of reactions selected from the core:", sum(fmc.find_reactions_from_pathway("Cholesterol metabolism")))