# Implements the WENO method used in Dumbser et al (DOI 10.1016/j.cma.2013.09.022)
# module Weno

# export coefficient_matrices, oscillation_indicator, weno


# include("array_manipulation.jl")
# include("Basis.jl")
# import Main.Basis: basis

function coefficient_matrices(N)
  # Generate linear systems governing the coefficients of the basis polynomials
  fHalfN = floor(N/2)
  cHalfN = ceil(N/2)
  M1 = zeros(N+1, N+1); M2 = zeros(N+1, N+1); M3 = zeros(N+1, N+1); M4 = zeros(N+1, N+1);
  ψ = basis(N)

  for p in 1:N+1, a in 1:N+1
    ψintp = integrate(ψ[p])
    M1[a,p] = ψintp(a-N)      - ψintp(a-N-1)
    M2[a,p] = ψintp(a)        - ψintp(a-1)
    M3[a,p] = ψintp(a-fHalfN) - ψintp(a-fHalfN-1)
    if isodd(N)
      M4[a,p] = ψintp(a-cHalfN) - ψintp(a-cHalfN-1)
    end
  end
  return factorize(M1), factorize(M2), factorize(M3), factorize(M4)
end

function oscillation_indicator(N)
  # Generate the oscillation indicator matrix
  Σ = zeros(N+1, N+1)
  ψ = basis(N)

  for a in 1:N
    ψdera = [derivative(ψp, a) for ψp in ψ]
    for p in 1:N+1, m in 1:N+1
      antiderivative = integrate(ψdera[p] * ψdera[m])
      Σ[p,m] += antiderivative(1) - antiderivative(0)
    end
  end
  return Σ
end

@inline function indicator(chΣT, w, n)
  # The indicator is w' * Σ * w, ie |w' * cholesky(Σ)'|^2
  # NB The Cholesky decomposition always exists as Σ is SPD
  # U = chΣT.U
  ret = 0.0
  for j in 1:n
    tmp = 0.0
    for i in j:n
      @inbounds tmp += w[i] * chΣT.U[i,j]
    end
    ret += tmp^2
  end
  return ret
end

function coeffs3(ret, N, w1,w2,w3, M1,M2,M3, λs,λc, chΣT, ε, r)
  # Calculate coefficients of basis polynomials and weights
  # A_ldiv_B!(M1, w1)
  # A_ldiv_B!(M2, w2)
  # A_ldiv_B!(M3, w3)
  w1 = M1 \ w1
  w2 = M2 \ w2
  w3 = M3 \ w3
  o1 = λs / (indicator(chΣT, w1, N+1) + ε)^r
  o2 = λs / (indicator(chΣT, w2, N+1) + ε)^r
  o3 = λc / (indicator(chΣT, w3, N+1) + ε)^r
  oSum = o1 + o2 + o3
  for i in 1:N+1
    @inbounds ret[i] = (o1*w1[i] + o2*w2[i] + o3*w3[i]) / oSum
  end
end

function coeffs4(ret, N, w1,w2,w3,w4, M1,M2,M3,M4, λs,λc, chΣT, ε, r)
  # Calculate coefficients of basis polynomials and weights
  # A_ldiv_B!(M1, w1)
  # A_ldiv_B!(M2, w2)
  # A_ldiv_B!(M3, w3)
  # A_ldiv_B!(M4, w4)
  w1 = M1 \ w1
  w2 = M2 \ w2
  w3 = M3 \ w3
  w4 = M4 \ w4
  o1 = λs / (indicator(chΣT, w1, N+1) + ε)^r
  o2 = λs / (indicator(chΣT, w2, N+1) + ε)^r
  o3 = λc / (indicator(chΣT, w3, N+1) + ε)^r
  o4 = λc / (indicator(chΣT, w4, N+1) + ε)^r
  oSum = o1 + o2 + o3 + o4
  for i in 1:N+1
    @inbounds ret[i] = (o1*w1[i] + o2*w2[i] + o3*w3[i] + o4*w4[i]) / oSum
  end
end

function weno(u, N, M1=nothing, M2=nothing, M3=nothing, M4=nothing, chΣT=nothing,
              λc=1e5, λs=1, r=8, ε=1e-14)
  # Find reconstruction coefficients of u to order N+1

  if M1==nothing
    M1,M2,M3,M4 = coefficient_matrices(N)
  end
  if chΣT==nothing
    Σ = oscillation_indicator(N)
    chΣT = cholesky(Σ)'
  end

  nx, ny, nz, nvar = size(u)
  fHalfN = floor(Int, N/2)
  cHalfN = ceil(Int, N/2)
  oddN = isodd(N)

  w1 = zeros(Float64, N+1)
  w2 = zeros(Float64, N+1)
  w3 = zeros(Float64, N+1)
  w4 = zeros(Float64, N+1)

  Wx = zeros(nx, ny, nz, N+1, nvar)
  Wx0 = add_ghost_cells(u, N, nx, ny, nz, nvar, 1)

  for v in 1:nvar, k in 1:nz, j in 1:ny
    w0 = view(Wx0, :, j, k, v)
    for i in 1:nx
      for ind  in  0:N
        @inbounds w1[ind+1] = w0[i+ind]
        @inbounds w2[ind+1] = w0[i+N+ind]
        @inbounds w3[ind+1] = w0[i+cHalfN+ind]
      end
      if oddN
        for ind  in  0:N
          @inbounds w4[ind+1] = w0[i+fHalfN+ind]
        end
        coeffs4(view(Wx,i,j,k,:,v), N, w1,w2,w3,w4, M1,M2,M3,M4, λs,λc, chΣT,ε,r)
      else
        coeffs3(view(Wx,i,j,k,:,v), N, w1,w2,w3, M1,M2,M3, λs,λc, chΣT,ε,r)
      end
    end
  end
  if ny==1 && nz==1
    return Wx
  end

  Wxy = zeros(nx, ny, nz, N+1, N+1, nvar)
  Wxy0 = add_ghost_cells(Wx, N, nx, ny, nz, nvar, 2)

  for v in 1:nvar, a in 1:N+1, k in 1:nz, i in 1:nx
    w0 = view(Wxy0, i, :, k, a, v)
    for j in 1:ny
      for ind  in  0:N
        @inbounds w1[ind+1] = w0[j+ind]
        @inbounds w2[ind+1] = w0[j+N+ind]
        @inbounds w3[ind+1] = w0[j+cHalfN+ind]
      end
      if oddN
        for ind  in  0:N
          @inbounds w4[ind+1] = w0[j+fHalfN+ind]
        end
        coeffs4(view(Wxy,i,j,k,a,:,v), N, w1,w2,w3,w4, M1,M2,M3,M4, λs,λc, chΣT,ε,r)
      else
        coeffs3(view(Wxy,i,j,k,a,:,v), N, w1,w2,w3, M1,M2,M3, λs,λc, chΣT,ε,r)
      end
    end
  end
  if nz==1
    return Wxy
  end

  Wxyz = zeros(nx, ny, nz, N+1, N+1, N+1, nvar)
  Wxyz0 = add_ghost_cells(Wxy, N, nx, ny, nz, nvar, 3)

  for v in 1:nvar, b in 1:N+1, a in 1:N+1, j in 1:ny, i in 1:nx
    w0 = view(Wxyz0, i, j, :, a, b, v)
    for k in 1:nz
      for ind  in  0:N
        @inbounds w1[ind+1] = w0[k+ind]
        @inbounds w2[ind+1] = w0[k+N+ind]
        @inbounds w3[ind+1] = w0[k+cHalfN+ind]
      end
      if oddN
        for ind  in  0:N
          @inbounds w4[ind+1] = w0[k+fHalfN+ind]
        end
        coeffs4(view(Wxyz,i,j,k,a,b,:,v), N, w1,w2,w3,w4, M1,M2,M3,M4, λs,λc, chΣT,ε,r)
      else
        coeffs3(view(Wxyz,i,j,k,a,b,:,v), N, w1,w2,w3, M1,M2,M3, λs,λc, chΣT,ε,r)
      end
    end
  end
  return Wxyz

end

"""
Add ghost cells which match boundary cells. This is equivalent to a homogeneous Neumann condition, I believe.
"""
function add_ghost_cells(arr, N, nx, ny, nz, nvar, d)
  # Extends the input array by N cells in the dth axis
  if d==1
    ret = zeros(nx+2N, ny, nz, nvar)
    ret[N+1:nx+N, :, :, :] = arr
    in0 = arr[1, :, :, :]
    in1 = arr[nx, :, :, :]
    for i = 1:N
      ret[i, :, :, :] = in0
      ret[nx+N+i, :, :, :] = in1
    end
    return ret

  elseif d==2
    ret = zeros(nx, ny+2N, nz, N+1, nvar)
    ret[:, N+1:ny+N, :, :, :] = arr
    in0 = arr[:, 1, :, :, :]
    in1 = arr[:, ny, :, :, :]
    # for i = 1:N
    #   ret[:, i, :, :] = in0
    #   ret[:, ny+N+i, :, :] = in1
    # end
    ret[:, 1:N, :, :, :] = in0
    ret[:, ny+N+1:ny+2N, :, :, :] = in1
    return ret

  else
    ret = zeros(nx, ny, nz+2N, N+1, N+1, nvar)
    ret[:, :, N+1:nz+N, :, :, :] = arr
    in0 = arr[:, :, 1, :, :]
    in1 = arr[:, :, nz, :, :, :]
    for i = 1:N
      ret[:, :, i, :, :, :] = in0
      ret[:, :, nz+N+i, :, :, :] = in1
    end
    return ret
  end
end

# end   # module Weno
