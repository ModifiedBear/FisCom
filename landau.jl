# Landau ginzburg hamiltonian + langevin dynamics

using DifferentialEquations
using Combinatorics: permutations, levicivita
using ImageFiltering: centered, imfilter
using GLMakie; GLMakie.activate!(float=true)
using MakieExtra
using ProgressBars
using Distributions

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

function diff_eq(u , p, t)
  # dp/dt = mu*F+eta
  μ,a,b,K,del_sq,dx = p
  return -μ*(a * u + 4*b*abs2.(u) .* u - K * imfilter(u, del_sq, "circular") / dx^2) # periodic boudnary conditions
  # return du

end

function diff_eq(u::Array{ComplexF64, 3} , p::Tuple, t::Float64)
  # dp/dt = mu*F+eta
  μ,a,b,K,del_sq,dx = p
  # return -μ*(a * u + 4*b*abs2.(u) .* u - K * imfilter(u, del_sq, "circular") / dx^2) # periodic boudnary conditions
  return -μ*(a * u + 4*b*abs2.(u) .* u - K * imfilter(u, del_sq, "circular") / dx^2) + 0.2*randn(ComplexF64, size(u)) # with noise
end

begin
  L          = 3.0
  mu         = 1.0
  a          = -1.0
  b          = 1.0 + 1.5im
  K          = 1.0



  Ndims      = (60,60)#, 20)#,51);
  Coords     = LinRange.(-L,L,Ndims);

  space_step = 1.0; # spatial width
  t0    = 0.0
  t1    = 5
  

  R = rand((-1,1),Ndims...); # initial field
  # R = [(z<0) for x in Coords[1], y in Coords[2], z in Coords[3]]
  # R = [exp(-(x^2+y^2+z^2) * 20) for x in Coords[1], y in Coords[2], z in Coords[3]] # initial field

  # θ = rand(Float64, Ndims...); # initial field
  V = copy(ComplexF64.(R))

  kernel_order = 2 # O(h^4)
  del = ComplexF64.(centered(get_laplace_kernel(length(Ndims),kernel_order))) # grad is not really complex, but for imfilter purposes it needs to

  params = (mu,a,b,K,del,space_step)


  MAXITER = 2500
  GC.gc()
  wave = copy(V)
  time_step = 0.01
  field = zeros(ComplexF64, MAXITER, Ndims...); # store wave
  # good old Euler method
  # for _i in ProgressBar(1:MAXITER)
  #   # wave = wave + time_step * (imfilter(wave, laplacian, "circular")/space_step^2 - abs2.(wave) .* wave + wave) # convolution
  #   wave = wave + diff_eq(wave, params, 0.0) * time_step
  #   field[_i,:,:,:] = wave # store wave
  # end



  prob = ODEProblem(diff_eq, V, (t0, t1), params)

  sol = solve(prob, Tsit5(), dt = 0.01)
  
  println("Interpolating...")
  time = LinRange(t0,t1, 1000)
  field = sol(time)
 
end



begin

  GC.gc()
  fig = Figure(size= (600,600))
  # sl = Slider(fig[0,1], range=range(1,MAXITER,step=1))
  ax = Axis3(fig[1,1]); ax.aspect=Ndims; 

  # psi = eachslice(field, dims=1)

  # vol = volume!(ax,real(field[1,:,:,:]), algorithm=:absorption, absorption=4.5f0)
  volre = volume!(ax, real(field[end,:,:,:]), algorithm=:iso, isovalue = 0.2, colormap=[:yellow, :orange], label = "Real")
  # volre = volume!(ax, abs2.(field[end,:,:,:]), algorithm=:absorption, absorption = 4, transparency=false, colormap=:cividis, label = "Real")
  # volim = volume!(ax, imag(field[end,:,:,:]), algorithm=:absorption, absorption = 4, transparency=true, colormap=:blues, label = "Real")
  volim = volume!(ax, imag(field[end,:,:,:]), algorithm=:iso, isovalue = 0.2, colormap=[:cyan, :purple], label = "Imag")
  # ax = Axis3(fig[1,1]); ax.aspect=(1,1,1); hm=volume!(ax,Coords..., wave_array[1], algorithm=:absorption, absorption=4.5f0, colormap=:balance)
  # cb = Colorbar(fig[2,1], hm, vertical=false, flipaxis=false, label=L"|\psi|^2"); cb.labelsize=24
  
  # ax.xticks = [
  Legend(fig[1,2], [PolyElement(color = :orange), PolyElement(color = :purple)], ["Real", "Imag"])
  hidedecorations!(ax)
  fig

  # Makie.lift(sl.value) do _i
  # _i = 2
  datare = eachslice(real.(field), dims=1)
  dataim = eachslice(imag.(field), dims=1)
  # while true
  record(fig, "landau_3D_noise_2.mp4", eachindex(datare), framerate=120) do _i

    
    volre[4][] = datare[_i]
    volim[4][] = dataim[_i]
    # _i += 5
    # ax.azimuth = _i / MAXITER * 4 * π / 10
    # sleep(0.01)
    # display(fig)

  end
  fig
end



# visualize
begin
  GC.gc()
  fig = Figure()

  ax = Axis(fig[1,1])
  data = real.(field)
  # data = [fill(zeros(size(V)),10); data...] # add some time before

  hm=heatmap!(ax,data[1],interpolate=false, colormap=:Spectral); ax.aspect=DataAspect()
  # vol = volume!(ax, real(field[end]), transparency=false, colormap=:Spectral)
  # # ax.aspect=(1,1,1)
  # ax.perspectiveness = 0.8
  # ax.elevation = pi/4 * 0.8
  # vol = contour!(ax, real(field[1]), transparency=false, colormap=:Spectral)
  # vol = volume!(ax, real(field[1]), algorithm = :absorption, absorption=4f0, colormap=:balance, transparency=false)

 

  cb = Colorbar(fig[1,2], hm, label="Re[m]")
  # cb.tickformat=EngTicks(:symbol)

  record(fig, "landau_2D_noise_2.mp4", eachindex(data), framerate=60) do _i
  # lift(sl.value) do _i
    hm[3][] = data[_i]
    # vol[4][] = data[_i]
    # ax.azimuth = _i / MAXITER 
    # _i += 100
    # sleep(0.00000000001)
  end
  fig
end
# heatmap(sol.u[10], colormap=:bam)