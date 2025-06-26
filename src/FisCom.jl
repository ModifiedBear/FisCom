module FisCom

# export oh2, oh4, oh6
include("./waves_FDTD.jl")
include("./Ising.jl")
using Permutations

# implement spatial gradient as convolution operator
function finite_diff_coefficient(_ord::Int64)
  if _ord == 2
    return [1,-2,1]
  elseif _ord == 4
    return [-1/12, 4/3, -5/2, 4/3, -1/12]
  elseif _ord == 6
    return [1/90, -3/20, 3/2, -49/18,3/2,-3/20,1/90]
  elseif _ord == 8
    return [-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]
  end
end

function get_laplace_kernel(_dim::Int64,_ord::Int64)
  base_ker = finite_diff_coefficient(_ord)
  kernel = zeros(fill(length(base_ker),_dim)...)
  mid_index = _ord÷2+1
  kernel[fill(mid_index,ndims(kernel)-1)...,:] .= base_ker # middle of middles
  
  perm_ind = collect(permutations(1:_dim))#[1:_dim-1:end]
  if _dim == 3
  # filter!(_i -> levicivita(_i) < 0, perm_ind)
  filter!(_i -> levicivita(_i) < 0, perm_ind)
  end
  kernel = mapreduce(_i -> permutedims(kernel,_i),+,perm_ind)
  

  #two_dim_ker += permutedims(two_dim_ker)
  return centered(kernel) # return centered version of kernel (0 = middle)
end

end