
module Weno

using FastGaussQuadrature
using Polynomials
using LinearAlgebra
using Polynomials

export coefficient_matrices, oscillation_indicator, weno

include("basis.jl")
include("weno_coeffs.jl")
end
