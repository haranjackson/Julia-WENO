# Weno.jl
An optimized Julia implementation of the WENO reconstruction algorithm of any order. Based on the work of Dumbser, Hidalgo, and Zanotti in *High Order Space-Time Adaptive WENO Finite Volume Schemes for Non-Conservative Hyperbolic Systems* (DOI 10.1016/j.cma.2013.09.022).

Implemented by Hari Jackson (@haranjackson) initially, updated for Julia 1.7 and reformulated as a package by Isaac Wheeler (@ickaser).

## Mathematical details

This implementation adds ghost cells *identical to the boundary cells*. This may not be appropriate for all problems.

## Usage
Your input data must take the form of an array `u` with shape `(nx,ny,nz,nvar)`.
`nvar` is the number of variables contained in each cell and `nx,ny,nz` are the number of cells in the *x,y,z*-axes, respectively.
If your grid is 1D or 2D, set `ny=1` and/or `nz=1`, as necessary.
Choose integer `N` where `N+1` is the desired order of accuracy.
The coefficients of the WENO reconstruction are obtained by calling:

```julia
julia> weno(u,N);
```

This returns an array with shape `(nx, ny, nz, N+1, N+1, N+1, nvar)`, where the first three dimensions correspond to physical space and the last corresponds to variables. Other dimensions are cell-by-cell coefficients for a nodal basis $\psi_l$ as described in the above paper. This polynomial basis is returned from `basis(N)`.

This call will construct the WENO coefficient matrices `M1,M2,M3,M4` and the oscillation indicator matrix `Σ` each time.
If many WENO reconstructions are required, it will be more computationally efficient to precalculate these entities:

```julia
julia> using LinearAlgebra;
julia> M1,M2,M3,M4 = coefficient_matrices(N);
julia> Σ = oscillation_indicator(N);
julia> chΣT = cholesky(Σ)';
julia> weno(u,N,M1,M2,M3,M4,chΣT);
```

`test/test.jl` contains a 1D example, `sin_cos_test`, which plots coefficients directly in space (since they correspond to interpolated points on the reconstruction).
