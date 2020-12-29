# MuFiM
MultiFrontal Factorization for General Sparse Matrix

## Installation


### MuFiM Installation

```
git clone --recursive https://github.com/YingzhouLi/MuFiM.git
cd MuFiM
matlab -nojvm -r "make;quit"
```

### Suggested Matlab toolbox

MuFiM could adopt both [meshpart] and Metis to bipartition the graph of a sparse matrix.
[meshpart] is a pure Matlab graph partitioning implementation, which is slower than METIS.
Therefore, it is suggested to install MetisMex toolbox before installing MuFiM.
The installation of MetisMex follows,

```
git clone https://github.com/YingzhouLi/metismex.git
cd metismex
matlab -nojvm -r "make;quit"
```

If metis package is not preinstalled, please find the detailed installation instructions at https://github.com/YingzhouLi/metismex.

[meshpart]: https://github.com/YingzhouLi/meshpart
