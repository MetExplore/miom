# Flux Balance Analysis (LP)

A simple Flux Balance Analysis can be easily defined and solved with any of the commercial or open-source solvers available:

```python
import miom
# Load  the iMM1865 mouse model from the repository
network = miom.mio.load_gem("@iMM1865")
target = "BIOMASS_reaction"
flux = (miom
        .load(network)
        .steady_state()
        .set_rxn_objective(target)
        .solve(verbosity=1)
        .get_fluxes(target))

print("Optimal flux is", flux)
```