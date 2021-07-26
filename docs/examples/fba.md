# Flux Balance Analysis (LP)

A simple Flux Balance Analysis can be easily defined and solved with any of the commercial or open-source solvers available:

```python
from miom import miom, load_gem
network = load_gem("https://github.com/pablormier/miom-gems/raw/main/gems/mus_musculus_iMM1865.miom")
target = "BIOMASS_reaction"
flux = (miom(network)
        .steady_state()
        .set_rxn_objective(target)
        .solve(verbosity=1)
        .get_fluxes(target))

print("Optimal flux is", flux)
```