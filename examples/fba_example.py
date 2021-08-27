import miom

# Load a model
m = miom.load_gem('@iMM1865')

target = 'BIOMASS_reaction'
opt_flux = (miom.load(m, solver=miom.Solvers.GLPK)
                .steady_state()
                .set_rxn_objective(target)
                .solve(verbosity=1)
                .get_fluxes(target))
        
# Optimal flux value
print(f"Optimal flux for {target} = {opt_flux}")
