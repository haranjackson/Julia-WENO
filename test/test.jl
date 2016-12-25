import PyPlot: plot

import Basis: nodes
import Weno: coefficient_matrices, oscillation_indicator, weno

function plot_weno(wh, var, N)
  n = size(wh)[1]
  x = float(repeat(0:n-1, inner=N+1))
  for i = 1:n
    x[1+(i-1)*(N+1) : i*(N+1)] += nodes(N)
  end
  y = zeros(n*(N+1))
  for i = 1:n
    for j = 1:N+1
      y[(i-1)*(N+1)+j] = wh[i, 1, 1, j, var]
    end
  end
  plot(x,y)
end

function main(u, N)
  M1,M2,M3,M4 = coefficient_matrices(N)
  Σ = oscillation_indicator(N)
  chΣT = chol(Σ)'
  @time weno(u, N, M1, M2, M3, M4, chΣT)
end

function sin_cos_test()
  u = rand(91,1,1,18)
  u[:,1,1,1] = sin(0:0.1:9)
  u[:,1,1,2] = cos(0:0.1:9)
  main(u,2)
  wh = main(u,2)
  plot_weno(wh,1,2)
  plot_weno(wh,2,2)
end
