using GLMakie; GLMakie.activate!(float=true) 
using ProgressBars
# using FisCom: FDTD

MAX_MEMORY = 5e9

# 2D example1

function oh6(vₙ, vₙ₊₁, α, κ, x, y,n)
  nx = length(x)
  ny = length(y)
  # better than the vectorized format
  bndry = 3;

  ffield = [exp(-5 * ((x)^2 + y^2)) for x in x, y in y]
  for ii in 1:n
    # you NEED to use copy(), this isn't python
    ω = 0.1
    #u_new[nx-2:nx-1,:] .= 10 .* sin(ω * ii)
    #M[ii,:,:] = copy(vₙ₊₁)
    #vₙ₊₁[nx-bndry-1:nx-bndry,div(ny,3):ny-div(ny,3)] .= 10 .* sin(ω * ii);
    vₙ₊₁ .= 10 .* cos(ω * ii) .* ffield;
    vₙ₋₁ = copy(vₙ);
    vₙ = copy(vₙ₊₁);
    for ii in 4:nx-3
      for jj in 4:ny-3
        vₙ₊₁[ii, jj] =α[ii,jj] * (  2  * vₙ[ii,jj-3]
                                  - 27 * vₙ[ii,jj-2]
                                  + 270* vₙ[ii,jj-1]
  
                                  + 2  * vₙ[ii-3,jj]
                                  - 27 * vₙ[ii-2,jj]
                                  + 270* vₙ[ii-1,jj]
                                  - 980* vₙ[ii,  jj]
                                  + 270* vₙ[ii+1,jj]
                                  - 27 * vₙ[ii+2,jj]
                                  + 2  * vₙ[ii+2,jj]
  
                                  + 270* vₙ[ii,jj+1]
                                  - 27 * vₙ[ii,jj+2]
                                  + 2  * vₙ[ii,jj+2]
                                  ) / 180 + 2 * vₙ[ii, jj] - vₙ₋₁[ii, jj];
        for kk in 1:bndry
          vₙ₊₁[kk, jj]     = vₙ[kk+1,   jj]   + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[kk+1,     jj]-vₙ[kk,     jj]); # x = 0 
          vₙ₊₁[nx-3+kk,jj] = vₙ[nx-4+kk , jj] + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[nx-4+kk,  jj]-vₙ[nx-3+kk,jj]); # x = N
          vₙ₊₁[ii, kk]     = vₙ[ii,  kk+1]    + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,     kk+1]-vₙ[ii,     kk]); # y = 0
          vₙ₊₁[ii,ny-3+kk] = vₙ[ii, ny-4+kk]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,  ny-4+kk]-vₙ[ii,ny-3+kk]);# y = N
        end
      end
    end
  end
  
  return vₙ₊₁
end  

begin
L = 5.0
nx,ny   = 91,91;
X,Y      = LinRange(-L,L,nx),LinRange(-L,L,ny);
speed      = 0.3;
space_step = 2.;# spatial width
time_step  = 0.5; # time step width (remember, dt < dx^2/2)?
velocities = ones(nx,ny) * speed # velocities 
KAPPA      = copy(velocities) * time_step/space_step  # 
ALPHA      = (velocities*time_step/space_step).^2; # wave equation coefficient 

V = zeros(nx, ny); # new
# p = Point2f.(L/2, -L/2)
# V = [exp(-2*((x+p[1])^2+(y+p[2])^2)) for x in X, y in Y]; # new
# V .= [exp(-2(x-2L/3)^2 + -0.1y^2) for x in X, y in Y]
W = zeros(nx, ny); # old


SIM_STEPS = 1400;
# wave_oh2 = FDTD.oh2(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=false); GC.gc()
# wave_oh4 = FDTD.oh4(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=false); GC.gc()
wave = oh6(V,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS); GC.gc()
end


begin
  GC.gc()
#  sl = SliderGrid(fig[2,1],(label="t", range=2:SIM_STEPS,startvalue=2))
#  hm = lift(sl.sliders[1].value) do _i
  fig = Figure(size=(500,400))
  ax = Axis(fig[1,1])
  hm = heatmap!(ax,X,Y, wave[2,:,:], colormap=:balance,interpolate=false, colorrange=(-1,1))
  ax.aspect = DataAspect()
  # lines!(ax, r_out * cos.(LinRange(0,2π,100)) .+ L/3, r_out * sin.(LinRange(0,2π,100)).+ L/3, color=:white, linewidth=1.5)
  # lines!(ax, r_out * cos.(LinRange(0,2π,100)) .- L/3, r_out * sin.(LinRange(0,2π,100)).- L/3, color=:white, linewidth=1.5)
  # vol_data = contour!(ax,X,Y,Z,abs.(volumes[2]), colormap=:turbo, levels=0.1:0.1:0.6)
#  vol_data = contour!(ax,X,Y,Z,volumes[2], transparency=false, colormap=:turbo, levels=0.1:0.1:0.2)
  for _i in 3:SIM_STEPS
  #  text!(ax, 0.0, 0.0, text="$_i")
    hm[3]=abs.(wave[_i,:,:])
    sleep(0.001)
    println(_i)
    display(fig)
  end
  #volume!(ax, X,Y,Z,volumes[_i])
#  end
  fig
end
