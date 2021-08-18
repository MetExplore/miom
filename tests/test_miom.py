from miom.miom import (
    Comparator, 
    ExtractionMode, 
    BaseModel,
    PythonMipModel,
    PicosModel
)
import pytest
from miom.mio import load_gem
import pathlib
import numpy as np


@pytest.fixture()
def gem():
    file = pathlib.Path(__file__).parent.joinpath("models", "example_r13m10.miom")
    return load_gem(str(file))

@pytest.fixture(params=[PythonMipModel, PicosModel])
def model(request, gem):
    # Get the different implementations and instantiate
    # them providing the test gem
    impl = request.param
    solver = 'cbc' if impl == PythonMipModel else 'glpk'
    return impl(miom_network=gem, solver_name=solver)

def prepare_fba(model, rxn=None, direction='max'):
    m = model.steady_state()
    if rxn is not None:
        m = m.set_rxn_objective(rxn, direction=direction)
    return m

def test_load_gem(model):
    assert model.network.num_reactions == 13


def test_fba_max(model):
    rxn = 'EX_i'
    f = (
        prepare_fba(model, rxn, direction='max')
        .solve()
        .get_fluxes(rxn)
    )
    assert np.isclose(f, 40/3)

def test_fba_min(model):
    rxn = 'EX_i'
    f = (
        prepare_fba(model, rxn, direction='min')
        .solve()
        .get_fluxes(rxn)
    )
    assert np.isclose(f, -20.0)


def test_subset_selection(model):
    m = prepare_fba(model, 'EX_h', direction='max')
    V, X = (
        m
        .solve()
        .set_fluxes_for('EX_h')
        .subset_selection(-1)
        .solve()
        .get_values()
    )
    assert np.sum(X) == model.network.num_reactions - 6

def test_network_selection_using_indicators(model):
    m = prepare_fba(model, 'EX_h', direction='max')
    network = (
        m
        .solve()
        .set_fluxes_for('EX_h')
        .subset_selection(-1)
        .solve()
        .select_subnetwork(
          mode=ExtractionMode.INDICATOR_VALUE,
          comparator=Comparator.LESS_OR_EQUAL,
          value=0.5
        )
        .network
    )
    assert network.num_reactions == 6

def test_network_selection_using_fluxes(model):
    m = prepare_fba(model, 'EX_h', direction='max')
    network = (
        m
        .solve()
        .set_fluxes_for('EX_h')
        .subset_selection(-1)
        .solve()
        .select_subnetwork(
          mode=ExtractionMode.ABSOLUTE_FLUX_VALUE,
          comparator=Comparator.GREATER_OR_EQUAL,
          value=1e-8
        )
        .network
    )
    assert network.num_reactions == 6


def test_mip_flux_consistency(model):
    V, X = (
        prepare_fba(model)
        .subset_selection(1)
        .solve()
        .get_values()
    )
    assert np.sum(X > 0.99) == model.network.num_reactions

def test_mip_flux_consistency_with_blocked_rxn(model):
    i, rxn = model.network.find_reaction("EX_a")
    rxn["lb"] = 0
    rxn["ub"] = 0
    V, X = (
        prepare_fba(model)
        .subset_selection(1)
        .solve()
        .get_values()
    )
    assert np.sum(X > 0.99) == model.network.num_reactions - 1

def test_subset_selection_custom_weights(model):
    # Large negative values
    c = [-100] * model.network.num_reactions
    # Positive weight only for R_f_i. Should not be selected
    # based on the objective function and the steady state constraints
    #i, _ = model.network.find_reaction('R_f_i')
    c[model.network.get_reaction_id('R_f_i')] = 1
    V, X = (
        prepare_fba(model)
        .subset_selection(c)
        .solve()
        .get_values()
    )
    # R_f_i has an indicator value of 1.0 only if the reaction is selected
    # since it's associated with a positive weight. The remaining reactions 
    # have an indicator value of 1.0 only if they are not selected (they are
    # associated with a negative weight). Since R_f_i should not be selected,
    # the expected number of ones in the indicator variables should be equal
    # to the number of reactions with negative weights that should not be
    # selected.
    assert np.sum(X > 0.99) == np.sum(np.array(c) < 0)

def test_activity_values(model):
    # Same problem as above, but now we use the activity values instead
    c = [-100] * model.network.num_reactions
    c[model.network.get_reaction_id('R_f_i')] = 1
    V, X = (
        prepare_fba(model)
        .subset_selection(c)
        .solve()
        .get_values()
    )
    active = sum(1 if abs(activity) == 1 else 0 for activity in model.variables.reaction_activity)
    inconsistent = sum(1 if activity != activity else 0 for activity in model.variables.reaction_activity)
    assert active == 0 and inconsistent == 0

def test_keep_rxn(model):
    # Large negative values
    c = [-100] * model.network.num_reactions
    # Positive weight only for R_f_i. Should not be selected
    # based on the objective function and the steady state constraints
    i = model.network.get_reaction_id('R_f_i')
    c[i] = 1
    V, X = (
        prepare_fba(model)
        .subset_selection(c)
        .keep('R_f_i') # Force to keep R_f_i anyway
        .solve()
        .get_values()
    )
    assert abs(V[i]) > 1e-8


    