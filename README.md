# Numerical computation of zigzag semigroup spectrum

This repository contains [Julia](https://julialang.org/) connected to the paper *Spectral analysis of the Zig-Zag process* by Joris Bierkens and Sjoerd Verduyn Lunel, 2019.

The code depends on the Julia packages QuadGK, PyPlot and SpecialFunctions, which may be installed from Julia using the command
```julia
using Pkg
Pkg.add(["PyPlot","QuadGK","SpecialFunctions"])
```
The spectrum for different target distributions can then be computed in the way displayed in the Jupyter Notebook file `compute_eigenvalues.ipynb`.
