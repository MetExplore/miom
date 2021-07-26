# iMAT: Integrative Metabolic Analysis Tool (MIP)

Example implementation of iMAT using MIOM. Note that this implementation supports also custom weights for the reactions. By default, weights for reactions in the "active" set are set to 1, and reactions in the "inactive" set are set to -1, so the objective function is exactly the same as the original iMAT.

!!! cite
    Shlomi, T., Cabili, M. N., Herrgård, M. J., Palsson, B. Ø., & Ruppin, E. (2008). Network-based prediction of human tissue-specific metabolism. [Nature biotechnology, 26(9), 1003-1010](https://www.nature.com/articles/nbt.1487).

```python
from miom import miom, Solvers
from miom.mio import load_gem
from miom.miom import Comparator, ExtractionMode

# Use the iHuman-GEM model
m = load_gem('https://github.com/pablormier/miom-gems/raw/main/gems/homo_sapiens_human1.miom')

# Add all the reactions from the Cholesterol pathway to the highly expressed set
RH = m.find_reactions_from_pathway("Cholesterol metabolism")
# Add reactions from the pyruvate metabolism to the lowly expressed set
RL = -1 * m.find_reactions_from_pathway("Pyruvate metabolism")
w = RH + RL
print("RH:", sum(RH), "RL:", sum(abs(RL)))

m = (miom(m, solver=Solvers.GUROBI)
     .setup(int_tol=1e-8, opt_tol=0.01, verbosity=1)
     .steady_state()
     .subset_selection(w)
     .solve(max_seconds=30)
     .select_subnetwork(
          mode=ExtractionMode.ABSOLUTE_FLUX_VALUE,
          comparator=Comparator.GREATER_OR_EQUAL,
          value=1e-8
     )
     .network)
print(m.num_reactions)

```