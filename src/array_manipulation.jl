function extend(arr, N, nx, ny, nz, nvar, d)
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
    in0 = arr[:, 1, :, :]
    in1 = arr[:, ny, :, :]
    for i = 1:N
      ret[:, i, :, :] = in0
      ret[:, ny+N+i, :, :] = in1
    end
    return ret

  else
    ret = zeros(nx, ny, nz+2N, N+1, N+1, nvar)
    ret[:, :, N+1:nz+N, :, :, :] = arr
    in0 = arr[:, :, 1, :]
    in1 = arr[:, :, nz, :]
    for i = 1:N
      ret[:, :, i, :] = in0
      ret[:, :, nz+N+i, :] = in1
    end
    return ret
  end
end
