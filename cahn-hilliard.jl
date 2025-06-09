using DifferentialEquations
using Combinatorics: permutations, levicivita
using ImageFiltering: centered, imfilter, Fill
using GLMakie; GLMakie.activate!(float=true)
using MakieExtra
using ProgressBars
using Distributions
using StatsBase: sample, Weights

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


function cahn_hilliard(u::Array{ComplexF64, 3} , p::Tuple, t::Float64)
  # dp/dt = mu*F+eta

  D, gamma, del_sq, dx, border = p

  # return -μ*(a * u + 4*b*abs2.(u) .* u - K * imfilter(u, del_sq, "circular") / dx^2) # periodic boudnary conditions
  mu = u.^3 - u - gamma * imfilter(u, del_sq, border) / dx^2

  return D * imfilter(mu, del_sq , border) / dx^2
end

function cahn_hilliard(u::Array{Float64, 3} , p::Tuple, t::Float64)
  # dp/dt = mu*F+eta

  D, gamma, del_sq, dx, border = p

  # return -μ*(a * u + 4*b*abs2.(u) .* u - K * imfilter(u, del_sq, "circular") / dx^2) # periodic boudnary conditions
  mu = u.^3 - u - gamma * imfilter(u, del_sq, border) / dx^2

  return D * imfilter(mu, del_sq , border) / dx^2 #+ 0.2*randn(Float64, size(u)) # with noise
end

begin
  D          = 1.0
  GAMMA      = 0.8
  L          = 3.0
  

  Ndims      = (60,60, 20)#, 20)#,51);
  Coords     = LinRange.(-L,L,Ndims);

  DX = 1.0; # spatial width
  # V  = rand((-1.0,1.0), Ndims...); # initial field
  # V = [sin(x)*cos(y) for x in Coords[1], y in Coords[2], z in Coords[3]] # initial field
  # V = rand((-1.0,1.0), Ndims...) # uniform mixing
  V = sample([-1.0,1.0], Weights([0.4, 0.6]), Ndims) # 60/40 mixing

  kernel_order = 6 # O(h^4)

  DEL = centered(get_laplace_kernel(length(Ndims),kernel_order))
  # del = ComplexF64.(centered(get_laplace_kernel(length(Ndims),kernel_order)))

  # pms = (D, GAMMA, DEL, DX, Fill(0.0, DEL))
  pms = (D, GAMMA, DEL, DX, "circular")

  MAXITER = 2500
  
  GC.gc()
  wave = copy(V)
  time_step = 0.005
  field = zeros(Float64, MAXITER, Ndims...); 

  # good old Euler method
  for _i in ProgressBar(1:MAXITER)
    
    wave = wave + cahn_hilliard(wave, pms, 0.0) * time_step
    if sum(isnan.(wave)) > 2
      println("NaN detected, stopping simulation")
      break
    end
    #mu = wave .^ 3 - wave - GAMMA * imfilter(wave, DEL, "circular") / DX^2
    #wave = wave + D * imfilter(mu, DEL , "circular") / DX^2 * time_step
    field[_i,:,:,:] = wave # store wave
  end


end

fig, ax,vol = volume(field[end,:,:,:], algorithm=:iso, transparency=true, isovalue=0.5)
volume!(ax, field[end,:,:,:], algorithm=:iso, transparency=true, isovalue=-0.5)
fig


mycat(a,b) = cat(a,b,dims=3)

bigvol = mycat(mycat(field[end, :,:,:],field[end,:,:,:]), field[end, :,:,:])

Makie.to_colormap(:balance)[end-10]
let
  GC.gc()

  data = eachslice(field,dims=1)

  fig = Figure(size= (600,600))
  # sl = Slider(fig[0,1], range=range(1,MAXITER,step=1))
  ax = Axis3(fig[1,1]); 

  vol_pos = volume!(ax, mycat(mycat(field[1,:,:,:],field[1,:,:,:]),field[1,:,:,:]), algorithm=:iso, isovalue=0.5, colormap=:balance)
  vol_neg = volume!(ax, mycat(mycat(field[1,:,:,:],field[1,:,:,:]),field[1,:,:,:]), algorithm=:iso, isovalue=-0.5, colormap=:balance)
  # volre = volume!(ax, abs2.(field[end,:,:,:]), algorithm=:absorption, absorption = 4, transparency=false, colormap=:cividis, label = "Real")
  # volim = volume!(ax, imag(field[end,:,:,:]), algorithm=:absorption, absorption = 4, transparency=true, colormap=:blues, label = "Real")
  # ax = Axis3(fig[1,1]); ax.aspect=(1,1,1); hm=volume!(ax,Coords..., wave_array[1], algorithm=:absorption, absorption=4.5f0, colormap=:balance)
  # cb = Colorbar(fig[2,1], hm, vertical=false, flipaxis=false, label=L"|\psi|^2"); cb.labelsize=24
  
  hidedecorations!(ax)

  Legend(fig[1,2], [PolyElement(color = Makie.to_colormap(:balance)[end-10]  ), 
                    PolyElement(color = Makie.to_colormap(:balance)[10])], ["Phase 1", "Phase 2"])

  # ax.xticks = [
  ax.elevation = 0.7
  fig

  # Makie.lift(sl.value) do _i
  # _i = 2
  # _i = 2
  # while true
  record(fig, "cahn_hilliard_unbalanced_PBC.mp4", eachindex(data), framerate=120) do _i

    
    vol_pos[4][] = mycat(mycat(field[mod(_i, MAXITER) + 1, :,:,:],field[mod(_i, MAXITER) + 1, :,:,:]),field[mod(_i, MAXITER) + 1, :,:,:])
    vol_neg[4][] = mycat(mycat(field[mod(_i, MAXITER) + 1, :,:,:],field[mod(_i, MAXITER) + 1, :,:,:]),field[mod(_i, MAXITER) + 1, :,:,:])

    _i += 4
    # ax.azimuth = _i / MAXITER * 4 * π / 10
    # sleep(0.01)
    # display(fig)

  end
end



volume(bigvol, algorithm=:iso, transparency=false, isovalue=0.5)
