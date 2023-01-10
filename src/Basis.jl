module Basis

export basis, nodes

import FastGaussQuadrature: gausslegendre
import Polynomials: Polynomial

function lagrange(x, w)
  """
  Taken from SciPy 0.18.1:

  Return a Lagrange interpolating polynomial.
  Given two 1-D arrays `x` and `w,` returns the Lagrange interpolating
  polynomial through the points ``(x, w)``.
  Warning: This implementation is numerically unstable. Do not expect to
  be able to use more than about 20 points even if they are chosen optimally.
  Parameters
  ----------
  x : array_like
      `x` represents the x-coordinates of a set of datapoints.
  w : array_like
      `w` represents the y-coordinates of a set of datapoints, i.e. f(`x`).
  Returns
  -------
      The Lagrange interpolating polynomial.
  """
  M = size(x, 1)[1]
  p = Polynomial([0])
  for j=1:M
    pt = Polynomial([w[j]])
    for k=1:M
      if k != j
        fac = x[j]-x[k]
        pt *= Polynomial([-x[k], 1.0]) / fac
      end
    end
    p += pt
  end
  return p
end

function nodes(N)
  # Returns Legendre-Gauss nodes, scaled to [0,1]
  return 0.5 .* (1 .+ gausslegendre(N+1)[1])
end

function basis(N)
  # Returns basis polynomials
  nodeArray = nodes(N)
  ret = Vector{Polynomial}()
  for i in 1:N+1
      basvec = zeros(N+1)
      basvec[i] = 1
      push!(ret, lagrange(nodeArray, basvec))
  end
  return ret
  # return [lagrange(nodeArray, eye(N+1)[:,i]) for i in 1:N+1]
end

end   # module Basis
