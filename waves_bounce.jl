using GLMakie; GLMakie.activate!(inline=false, float=true)
using FisCom: FDTD

# 2D example1
begin
L          = 5.0
nx,ny      = 91,91;
X,Y        = LinRange(-L,L,nx),LinRange(-L,L,ny);
speed      = 0.2;
space_step = 1.0;# spatial width
time_step  = 0.5; # time step width (remember, dt < dx^2/2)?
velocities = ones(nx,ny) * speed # velocities 
KAPPA      = copy(velocities) * time_step/space_step  # 
ALPHA      = (velocities*time_step/space_step).^2; # wave equation coefficient 

U = [exp(- 4(x)^2 -4(y)^2) for x in X, y in Y]
V = zeros(nx, ny); # new
# V = [exp(-5*((x-L/3)^2+(y-L/3)^2)) for x in X, y in Y]; # new
# V[:,4] .= [exp(-5*x^2) for x in X]

# V = 
W = zeros(nx, ny); # old

SIM_STEPS = 1200;
# wave_oh2 = FDTD.oh2(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
# wave_oh4 = FDTD.oh4(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
wave_oh6 = FDTD.oh6(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=false); GC.gc()
end

begin
  L          = 5.0
  nx         = 91;
  X          = LinRange(-L,L,nx);
  speed      = 0.2;
  space_step = 1.0;# spatial width
  time_step  = 0.5; # time step width (remember, dt < dx^2/2)?
  velocities = ones(nx) * speed # velocities 
  KAPPA      = copy(velocities) * time_step/space_step  # 
  ALPHA      = (velocities*time_step/space_step).^2; # wave equation coefficient 
  
  U = [exp(-x^2) for x in X]
  V = zeros(nx); # new
  # V = [exp(-5*((x-L/3)^2+(y-L/3)^2)) for x in X, y in Y]; # new
  # V[:,4] .= [exp(-5*x^2) for x in X]
  
  # V = 
  W = zeros(nx); # old
  
  SIM_STEPS = 4200;
  # wave_oh2 = FDTD.oh2(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
  # wave_oh4 = FDTD.oh4(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
  wave_oh6 = FDTD.oh6(U,copy(V),ALPHA, KAPPA, X, SIM_STEPS;_absorbing=true); GC.gc()
end

heatmap(wave_oh6)


begin
  GC.gc()
  MAX_VAL = maximum(wave_oh2)
#  sl = SliderGrid(fig[2,1],(label="t", range=2:SIM_STEPS,startvalue=2))
#  hm = lift(sl.sliders[1].value) do _i
  fig = Figure(size=(800,400))
  ax = [Axis(fig[1,1])]
  # plot_oh2 = heatmap!(ax[1],X,Y,abs.(wave_oh2[2,:,:]),colormap=:balance,interpolate=false,colorrange=(0,1))
  # plot_oh4 = heatmap!(ax[2],X,Y,abs.(wave_oh4[2,:,:]),colormap=:balance,interpolate=false,colorrange=(0,1))
  # plot_oh6 = heatmap!(ax[1],X,Y,abs2.(wave_oh6[2,:,:]),colormap=:dense,interpolate=false,colorrange=(0,10))
  plot_oh6 = lines!(ax[1], wave_oh6[2,:], color=:black)
  # vol_data = contour!(ax,X,Y,Z,abs.(volumes[2]), colormap=:turbo, levels=0.1:0.1:0.6)
#  vol_data = contour!(ax,X,Y,Z,volumes[2], transparency=false, colormap=:turbo, levels=0.1:0.1:0.2)
  [ax.aspect=DataAspect() for ax in ax]
  [hidedecorations!(ax) for ax in ax]
  ylims!(ax[1], -8,8)
  for _i in 3:SIM_STEPS
  #  text!(ax, 0.0, 0.0, text="$_i")
    # plot_oh2[3][]=abs.(wave_oh2[_i,:,:])
    # plot_oh4[3][]=abs.(wave_oh4[_i,:,:])

    plot_oh6[1]=wave_oh6[_i,:]
    sleep(0.001)
    println(_i)
    display(fig)
  end
  #volume!(ax, X,Y,Z,volumes[_i])
#  end
  fig
end
