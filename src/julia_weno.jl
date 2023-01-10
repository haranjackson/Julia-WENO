
module Weno

using FastGaussQuadrature
using Polynomials
using LinearAlgebra
import Polynomials: integrate, derivative

export coefficient_matrices, oscillation_indicator, weno
# module julia_weno
# include("add_ghost_cells.jl")
include("basis.jl")
include("weno.jl")
end
