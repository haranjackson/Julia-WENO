# Julia-WENO
An optimized Julia implementation of the WENO reconstruction algorithm, up to order 5, based on the work of Dumbser, Hidalgo, and Zanotti in *High Order Space-Time Adaptive WENO Finite Volume Schemes for Non-Conservative Hyperbolic Systems* (DOI 10.1016/j.cma.2013.09.022).

## Usage
Your input data must take the form of an array *u* with shape (*nx,ny,nz,nvar*) where *nvar* is the number of variables contained in each cell, and *nx,ny,nz* are the number of cells in the *x,y,z* axes, respectively. If your grid is 1D or 2D, set *nz=1* and *ny=1*, as necessary. Choose integer *N* where *N+1* is the desired order of accuracy, up to 5.

The function *main* contained in *weno.jl* takes *u* and *N* as inputs, and returns an array of WENO coefficients with shape (*nx,ny,nz,N+1,nvar*), (*nx,ny,nz,N+1,N+1,nvar*), or (*nx,ny,nz,N+1,N+1,N+1,nvar*), depending on whether your data is 1D, 2D, or 3D, respectively. The function *test* in *weno.jl* can be run to calculate and plot the WENO reconstruction of *sin(x)* on *1<x<10*.
