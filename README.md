<div align="center">
    <a href="https://metexplore.github.io/miom"><img align="center" src="https://github.com/MetExplore/miom/raw/development/docs/assets/img/miom_v1.png" alt="MIOM" title="MIOM: Mixed Integer Optimization for Metabolism"/>
    </a>
</div>

<div align="center">

[![Try It Online](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1JAOEHLlRCW8GziIpBqkFwJL2ha3OEOWJ?usp=sharing)
[![PyPI](https://img.shields.io/pypi/v/miom?label=PyPI&logo=pypi&logoColor=lightgrey)](https://pypi.org/project/miom/)
![Tests](https://github.com/metexplore/miom/actions/workflows/unit-tests.yml/badge.svg)
[![Downloads](https://pepy.tech/badge/miom)](https://pepy.tech/project/miom)

</div>

<h4 align="center">
<a href="https://metexplore.github.io/miom">https://metexplore.github.io/miom</a>
</h4>

---

__MIOM__ (Mixed Integer Optimization for Metabolism) is a python library for creating and solving complex optimization problems using genome-scale metabolic networks, in just a few lines. 

MIOM offers a high-level API that leverages the power of modern Mixed Integer Optimization (MIO) solvers to easily define steady-state metabolic optimization problems, from simple Flux Balance Analysis (FBA) simulations, to more complex problems, such as sparse FBA or context-specific reconstruction algorithms, and solve them the __required level of optimality__. It supports many free and commercial solvers thanks to the integration with the [PICOS](https://picos-api.gitlab.io/picos/) and the [Python-MIP](https://www.python-mip.com/). It is also compatible and complementary to [cobrapy](https://opencobra.github.io/cobrapy/).

Here is a simple example of the differences of implementing from scratch a simple Flux Balance Analysis (FBA), which is typically implemented as a Linear Program (LP) problem, and the Sparse FBA problem which is an Integer Programmming (IP) problem that solves the same FBA problem but minimizing the number of active reactions:

> NOTE: This library is functional but still in a very early stage. API is still not stable and might be changed in future versions.


```python
import miom

# Download the Recon3D metabolic network and find the maximum flux value
# through the biomass_reaction
model = (miom
        .load('@BiGG/Recon3D.miom')
        .steady_state()
        .set_rxn_objective('biomass_reaction')
        .solve(verbosity=1))

print("Optimal flux:", model.get_fluxes('biomass_reaction'))
print("Number of active reactions:", sum(abs(model.get_fluxes()) > 1e-8))

# Transform the previous FBA optimization problem into a Sparse FBA problem (MIP).
# Solving this problem returns the minimum number of active reactions preserving
# the optimal flux through the biomass_reaction
V, X = (model
        .setup(opt_tol=0.05)
        .set_fluxes_for('biomass_reaction')
        .subset_selection(-1)
        .solve(verbosity=1)
        .get_values())

print("Number of active reactions:", sum(abs(V) > 1e-8))
```

## Installation

### Minimal install

By default, MIOM comes with support for COIN-OR CBC solver and GLPK using the swiglpk wrapper. To install MIOM with minimal dependencies, run:

```
pip install miom
```

### Full install

The full install option of `miom` comes also with the interfaces for [Gurobi](https://www.gurobi.com/downloads) and [Mosek](https://www.mosek.com/downloads/). Also, in order to be able to import metabolic networks in different formats (SBML, Matlab, YAML, etc), you need to install the additional `cobra` and `scipy` packages:

```
pip install miom[all] cobra scipy
```

CPLEX is also supported, but requires a valid license. To install MIOM with CPLEX support, follow the instructions on the [CPLEX page](https://www.ibm.com/docs/en/icos/12.8.0.0?topic=cplex-setting-up-python-api) to install your current version of cplex in your python environment.


## Quick start

### Importing a Genome-Scale Metabolic Network (GEMs)

First step is to import a Genome-Scale Metabolic Network. The method `load_gem()` can be used to import a network from a local file or a URL. This method requires the [cobrapy](https://cobrapy.readthedocs.io/en/latest/) and [scipy](https://www.scipy.org/) packages to import networks in the SBML, YAML, JSON, and Matlab formats (see [Full install](#full-install)). Here is an example of importing the [BiGG Recon3D](http://bigg.ucsd.edu/models/Recon3D) network:

```python
import miom
network = miom.load_gem('https://github.com/SBRG/bigg_models_data/raw/master/models/Recon3D.mat')
print("Number of reactions in the network:", network.num_reactions)
```

By default, MIOM uses its own format (.miom) which is a lightweight, compressed and portable format based on numpy structured arrays to store only the essential information required to perform common tasks. To improve reproducibility of experiments, a repository of already converted models is available at [pablormier/miom-gems](https://github.com/pablormier/miom-gems/). Any model from this repository can be imported using shorthand links in the following forma: `@relative/path/to/model`. For example, in order to import the [Human-GEM v1.9.0](https://github.com/pablormier/miom-gems/tree/main/gems/SysBioChalmers/Human-Gem/v1.9.0), you only need to run:

```python
network = miom.load_gem('@SysBioChalmers/Human-Gem/v1.9.0')
```

This is a very convenient way of importing models and sharing **reproducible experiments**, making sure that other users are going to test the same models.

### Managing Metabolic Networks

MIOM is not intended to provide functions for creating and managing metabolic networks, since there are already other libraries for this purpose (e.g. [cobrapy](https://opencobra.github.io/cobrapy/)). If you rely on `cobrapy` for model manipulation, you can still use MIOM for the optimization after you prepare your metabolic network, as MIOM is nicely integrated with COBRA:

```python
from cobra.test import create_test_model
network = create_test_model("textbook")
medium = network.medium
medium["EX_o2_e"] = 0.0
network.medium = medium

# Use MIOM for optimization
flux = (miom
        .load(miom.mio.cobra_to_miom(network))
        .steady_state()
        .set_rxn_objective('Biomass_Ecoli_core')
        .solve()
        .get_fluxes('Biomass_Ecoli_core'))
```

Metabolic Networks in MIOM are represented by a `miom.mio.MiomNetwork` class. This class encodes the network into three matrices: `R` (reactions), `S` (stoichiometry), and `M` (metabolites). The `R` matrix is a numpy structured array that contains the list of the reactions defined in the network, including the reaction ID, the reaction name, the reaction bounds, the subsystem and the Gene-Protein-Rules:

```python
>>> miom.load_gem('@BiGG/Recon3D.miom').R[:5]

array([('10FTHF5GLUtl', '5-Glutamyl-10Fthf Transport, Lysosomal', 0., 1000., 'Transport, lysosomal', ''),
       ('10FTHF5GLUtm', '5-Glutamyl-10Fthf Transport, Mitochondrial', 0., 1000., 'Transport, mitochondrial', ''),
       ('10FTHF6GLUtl', '6-Glutamyl-10Fthf Transport, Lysosomal', 0., 1000., 'Transport, lysosomal', ''),
       ('10FTHF6GLUtm', '6-Glutamyl-10Fthf Transport, Mitochondrial', 0., 1000., 'Transport, mitochondrial', ''),
       ('10FTHF7GLUtl', '7-Glutamyl-10Fthf Transport, Lysosomal', 0., 1000., 'Transport, lysosomal', '')],
      dtype=[('id', 'O'), ('name', 'O'), ('lb', '<f8'), ('ub', '<f8'), ('subsystem', 'O'), ('gpr', 'O')])
```

Basic manipulation, for example, changing the bounds of a reaction, can be done just by modifying the corresponding attribute of the reaction:

```python
rxn = network.R[0]
rxn['lb'] = 0
rxn['ub'] = 10
```

Numpy structured arrays make it easy to modify properties of the network. For example, in order to change the lower bound of all reactions in the network at the same time, you only need to do:

```python
network.R['lb'] = -1000
```

Structured arrays also support indexing:

```python
import numpy as np
# List all reactions with a negative lower bound
network.R[network.R['lb'] < 0]
# Get the indexes of the reactions with a negative lower bound
idxs = np.flatnonzero(network.R['lb'] < 0)
```

### Constraint-based optimization problems

MIOM can be seen as a high-level API for creating constraint-based optimization problems by applying successive composable operations that add new constraints to the problem. The basic operation is the `steady_state()` method, which adds the system of equations `S * V = 0` to the optimization model, where `S` is the stoichiometric matrix and `V` is the vector of fluxes:

```python
model = (miom
        .load('@BiGG/Recon3D.miom')
        .steady_state())
```

Now the optimization model contains the system of equations for the steady state condition of the fluxes. Now, it can be extended by adding new constraints or changing properties of the model. For example, in order to find the maximum flux of the `biomass_reaction` (limiting the max production to 10) can be done as follows:

```python
flux = (model
        .set_rxn_objective('biomass_reaction')
        .set_flux_bounds('biomass_reaction', max_flux=10.0)
        .solve()
        .get_fluxes('biomass_reaction'))
```

More complex constraints can be added by using the `add_constraint()` method. This method accepts an affine expression involving the variables in the model. Here is an example of how to perform FBA again but with a constraint that ensures that the sum through all non-reversible reactions is less than or equal to 10:

```python
import numpy as np

# Load the Recon3D model from the MIOM repository
network = miom.load_gem('@BiGG/Recon3D.miom')
# The objective is to maximize the flux through the biomass reaction
model = miom.load(network).steady_state().set_rxn_objective('biomass_reaction')
# Find reactions that are not reversible (i.e. can have only positive fluxes)
non_reversible = np.flatnonzero(network.R['lb'] >= 0)
# Add the constraint that the sum through all non-reversible reactions is less than or equal to 10
# (note that both PICOS and Python-MIP are symbolic, and expressions involving the variables are allowed)
constraint = sum(model.variables.fluxvars[i] for i in non_reversible) <= 10
# Add to the model and solve 
model.add_constraint(constraint).solve(verbosity=1)
```
```
Starting solution of the Linear programming problem using Dual Simplex

Coin0506I Presolve 2242 (-3594) rows, 6144 (-4456) columns and 33659 (-12128) elements
Clp0014I Perturbing problem by 0.001% of 1.9797539 - largest nonzero change 1.4130091e-06 ( 7.1372964e-05%) - largest zero change 1.412981e-06
Clp0006I 0  Obj 0.0015868455 Primal inf 3303572.7 (1175) Dual inf 1.9797525 (1)
Clp0000I Optimal - objective value 9.062223
Coin0511I After Postsolve, objective 9.062223, infeasibilities - dual 0 (0), primal 0 (0)
Clp0032I Optimal objective 9.062223036 - 1675 iterations time 0.142, Presolve 0.07
```

### Mixed Integer Optimization

Many constraint-based metabolic problems consist of selecting a subset of reactions from a generic metabolic network, subject to some biological-based constraints. For example, Sparse FBA consists of selecting the minimum number of active reactions that lead to a desired optimal flux. Context-specific network reconstruction methods such as iMAT, Fastcore, mCADRE, MBA, etc., have also the same purpose: select a subset of reactions from a network based on some experimental data and some objective function. All these problems are instances of a more general problem, known as **best subset selection problem**, which can be modelled using *Mixed Integer Optimization (MIO)*. However, instances of this problem are usually hard to solve, and some implementations like Fastcore or MBA use different heuristics or relaxations to find sub-optimal solutions in short time. 

Unfortunately, there are a few problems with these type of approximations. First, we don't know how good is the solution obtained, and in practice many studies that use those methods assume that the solution they obtained is the best one. Second, even though most of these techniques solve the same type of problem, they look very different or hard to understand for practitioners. Third, many of these implementations do not exploit the potential of the modern MIO solvers which are nowadays able to solve very large problems in a short time, and to adjust the desired level of optimality, which gives the user an idea of how good is the solution obtained.

MIOM incorporates specific methods to easily model and solve such problems. Just by calling the method `subset_selection()`, which takes as input a weight for all reactions or a list of weights for each reaction, it transforms the current problem into a *best subset selection* problem, in which the objective function is the sum of the weights of the reactions in the solution (where positive weighted reactions contribute to the objective function if they are selected, and negative weighted reactions contribute to the objective function if they are not selected).

This simple method makes it possible to model a wide variety of complex constraint-based optimization problems. See for example how easy it is to implement the exact and generalized version of the Fastcore with MIOM:

```python
import miom
import numpy as np

# Use the flux-consistent subnetwork (fcm) of the Human1 GEM model
# NOTE: Fastcore requires that reactions in the core are flux consistent,
# otherwise the problem would be infeasible. 
m = miom.load_gem('@homo_sapiens_human1_fcm.miom')
# Select reactions from the cholesterol metabolism as the core reactions to keep
core_rxn = m.find_reactions_from_pathway("Cholesterol metabolism")
# Assign a negative weight for reactions not in the core
weights = -1 * np.ones(m.num_reactions)
weights[core_rxn == 1] = 1

fmc = (miom
        .load(m)
        .setup(opt_tol=0.05, verbosity=1)
        .steady_state()
        .subset_selection(weights)
        .keep(core_rxn == 1)
        .solve()
        .select_subnetwork(
            mode=miom.ExtractionMode.ABSOLUTE_FLUX_VALUE,
            comparator=miom.Comparator.GREATER_OR_EQUAL,
            value=1e-8
        )
        .network
)
print(fmc.num_reactions)
```

## Advantages

* __It's flexible:__ MIOM uses the [PICOS](https://picos-api.gitlab.io/picos/) and the [Python-MIP](https://www.python-mip.com/) libraries, which means you can use any solver supported by those libraries.
* __It's easy to extend:__ MIOM is written in pure python, so you can easily extend it to solve more complex optimization problems.
* __It makes the problem explicit:__ MIOM uses a declarative way to express the problem, so you can easily read and understand what you are solving and differences between algorithms.
* __It's fast:__ MIOM leverages the power of MIO solvers to solve complex optimization problems. You can control the quality and speed of the solutions for your problem and get better solutions that the approximations (LP) of the original problem available in other constraint-based modeling libraries.
* __It's lightweight:__ The library has a small number of dependencies, so it's easy to install and distribute also in HPC environments.
* __It includes a compressed GEM format__: MIOM can load and save the minimal information of the metabolic networks required for performing simulations into a compressed file compatible with numpy. The small size of the files allows you to quickly run online experiments so other people can reproduce your results. It also supports SBML and matlab formats if `cobratoolbox` is installed.
* __It's open-source:__ MIOM is open-source and free to use. You can contribute to the development of MIOM by forking the repository and sending pull requests.

## Documentation

The documentation of the library can be found at https://metexplore.github.io/miom/

## How to cite

Manuscript in progress

## License

GNU General Public License v3.0
