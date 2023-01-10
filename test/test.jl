import PyPlot: plot
using LinearAlgebra


function plot_weno(wh, var, N)
  n = size(wh, 1)
  x = repeat(0.0:n-1, inner=N+1)
  x += repeat(nodes(N), outer=n)
  y = reshape(wh[:, 1, 1, :, var]', :)
  plot(x,y)
end

function main(u, N)
  M1,M2,M3,M4 = coefficient_matrices(N)
  Σ = oscillation_indicator(N)
  chΣT = cholesky(Σ)'
  @time weno(u, N, M1, M2, M3, M4, chΣT)
end

function sin_cos_test(n=20)
  u = rand(n,1,1,2)
  xs = range(0, 9, length=n)
  u[:,1,1,1] = sin.(xs)
  u[:,1,1,2] = cos.(xs)
  main(u,2)
  wh = main(u,2)
  plot_weno(wh,1,2)
  plot_weno(wh,2,2)
end
