<div align="center">
    <a href="https://metexplore.github.io/miom"><img align="center" src="https://github.com/MetExplore/miom/raw/development/docs/assets/img/miom_v1.png" alt="MIOM" width="545" title="MIOM: Mixed Integer Optimization for Metabolism"/>
    </a>
</div>

<div align="center">

[![Try It Online](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1JAOEHLlRCW8GziIpBqkFwJL2ha3OEOWJ?usp=sharing)
[![PyPI](https://img.shields.io/pypi/v/probnum?style=flat-square&label=PyPI&logo=pypi&logoColor=white)](https://pypi.org/project/miom/)
![Tests](https://github.com/metexplore/miom/actions/workflows/unit-tests.yml/badge.svg)

</div>

---

https://metexplore.github.io/miom

__MIOM__ (Mixed Integer Optimization for Metabolism) is a python library for creating and solving complex optimization problems using genome-scale metabolic networks, in just a few lines. 

MIOM offers a high-level API that leverages the power of modern Mixed Integer Optimization (MIO) solvers to easily define steady-state metabolic optimization problems, from simple Flux Balance Analysis (FBA) simulations, to more complex problems, such as sparse FBA or context-specific reconstruction problems, and solve them the __required level of optimality__.

MIOM uses the [PICOS](https://picos-api.gitlab.io/picos/) and the [Python-MIP](https://www.python-mip.com/) libraries to build and solve the optimization problems using many commercial, academic and free solvers. It is also compatible and complementary to [cobrapy](https://opencobra.github.io/cobrapy/).

Here is an example of how to implement FBA and Sparse FBA to maximize flux through the biomass reaction in the Recon3D model with MIOM:

```python
import miom

# Download the Recon3D metabolic network and perform FBA
model = (miom
        .load('@BiGG/Recon3D.miom')
        .steady_state()
        .set_rxn_objective('biomass_reaction')
        .solve(verbosity=1))

print("Optimal flux:", model.get_fluxes('biomass_reaction'))
print("Number of active reactions:", sum(abs(model.get_fluxes()) > 1e-8))

# Transform the previous FBA optimization into a Sparse FBA optimization:
V, X = (model
        .setup(opt_tol=0.05)
        .set_fluxes_for('biomass_reaction')
        .subset_selection(-1)
        .solve(verbosity=1)
        .get_values())

print("Number of active reactions:", sum(abs(V) > 1e-8))
```

> NOTE: This library is functional but still in a very early stage. API is still not stable and might be changed in future versions.

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

The first step is to import a GEM (a genome-scale model) using the `miom.mio.load_gem` method. The method accepts a local file or an URL. SBML files and Matlab files are supported through `cobra` and `scipy` packages. Please make sure that these packages are correctly installed before trying to import models in these formats (see [Full install](#full-install)).

Once these packages are available, you can import a model using:

```python
import miom
network = miom.mio.load_gem("https://github.com/SysBioChalmers/Human-GEM/raw/main/model/Human-GEM.mat")
print("Number of reactions in the network", network.num_reactions)
```

MIOM includes its own lightweight and portable format (.miom) for loading and storing metabolic networks, which does not require any additional dependency. This models uses `numpy` compressed structured arrays to store the basic information of the metabolic networks required for performing simulations. This format is even smaller and more portable than the Matlab files. A repository of converted models is available at https://github.com/pablormier/miom-gems.

To make things even easier, the method `load_gem` can import any model from this repository using the relative path, for example:

```python
human1 = miom.mio.load_gem("@SysBioChalmers/Human-GEM.miom")
recon3d = miom.mio.load_gem("@BiGG/Recon3D.miom")
```

Here is an example of how to load a metabolic network and maximize the flux through a target reaction using FBA, and then how to modify the original problem to implement the sparse FBA problem adding only a few lines to the original problem:

```python
import miom

network = miom.mio.load_gem("@mus_musculus_iMM1865.miom")
target_rxn = "BIOMASS_reaction"

# Create the optimization problem with miom and solve
model = (miom
        .load(network)
        .steady_state()
        .set_rxn_objective(target_rxn)
        .solve(verbosity=1))

print("Optimal flux:", model.get_fluxes(target_rxn), "mmol/(h·gDW)")

# Show reactions with non-zero flux
V, X = model.get_values()
print("Number of reactions active reactions:", sum(abs(V) > 1e-8))
```

```
Optimal flux: 798.8110517749975 mmol/(h·gDW)
Number of active reactions: 2549
```

Now, modify the original problem to solve the sparse FBA problem, minimizing the number of reactions with non-zero flux that can lead to the optimal possible flux through the target reaction. This can be easily done by transforming the FBA problem into a subset selection problem, where each reaction has a negative weight and the goal is to remove as many negative weighted reactions as possible. Note that since each reaction has the same weight (-1), all reactions are equally important in the optimization problem:

> Note: To better understand the meaning of each step, please read the documentation of the [BaseModel class](https://metexplore.github.io/miom/references/miom/#miom.miom.BaseModel), and the complete example in [examples/sparse_fba.py](https://metexplore.github.io/miom/examples/sparse_fba).    
   

```python
V, X = (model
        # Set the MIP Gap tolerance to 5%, using the default solver
        # (COIN-OR CBC included with the Python-MIP lib)
        .setup(opt_tol=0.05, verbosity=1)
        # Set the fluxes to the optimal value found
        .set_fluxes_for('BIOMASS_reaction')
        # Transform to a subset selection problem
        # adding a negative weight to all reactions
        # (to remove them from the final solution)
        .subset_selection(-1)
        # Solve the MIO problem
        .solve()
        # Get continuos vars (fluxes) and binary vars
        .get_values())

# Show reactions with non-zero flux
print("Number of active reactions:", sum(abs(V) > 1e-8))
```

```
Number of active reactions: 404
```

Solving this problem with default COIN-OR CBC solver returns a solution with 404 active reactions (much less than the 2549 reactions obtained with FBA, and less than the 433 reactions returned by the CappedL1 approximation in the [sparseFBA](https://opencobra.github.io/cobratoolbox/stable/modules/analysis/sparseFBA/index.html) implementation in Matlab), with a relative gap between the lower and upper objective bound below 5% (as indicated in the setup method):

```
Cbc0011I Exiting as integer gap of 122.04538 less than 1e-10 or 5%
Cbc0001I Search completed - best objective -10208, took 0 iterations and 0 nodes (28.34 seconds)
Cbc0035I Maximum depth 0, 0 variables fixed on reduced cost
Total time (CPU seconds):       60.79   (Wallclock seconds):       60.79
```

The concise API provided by MIOM makes everything explicit: the sparse FBA problem can be implemented as a best subset selection problem of reactions (minimize the number of reactions with non-zero flux) subject to the steady-state constraint and the optimality constraint of the flux for the target reaction (in this case the `BIOMASS_reaction`). Using this formulation, you can take advantage of modern solvers like CPLEX, GUROBI, MOSEK, COIN-OR CBC (among others) to obtain an optimal or an approximate solution, controlled by the `opt_tol` parameter.

To use other solvers, you only need to provide the specific solver to the `miom` method, for example:

```python
model = (miom
        .load('@BiGG/Recon3D.miom', solver=miom.Solvers.GLPK)
        .steady_state()
        .set_rxn_objective('biomass_reaction')
        .solve(verbosity=1))
```

## Advantages

* __It's easy to use:__ MIOM uses the [PICOS](https://picos-api.gitlab.io/picos/) and the [Python-MIP](https://www.python-mip.com/) libraries, which means you can use any solver supported by those libraries.
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
