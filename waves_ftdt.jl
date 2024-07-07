using GLMakie
using FisCom: FDTD

# 2D example
begin
L = 5.0
nx,ny   = 91,91;
X,Y      = LinRange(-L,L,nx),LinRange(-L,L,ny);
speed      = 0.2;
space_step = 1.0;# spatial width
time_step  = 1.0; # time step width (remember, dt < dx^2/2)?
velocities = ones(nx,ny) * speed # velocities 
KAPPA      = copy(velocities) * time_step/space_step  # 
ALPHA      = (velocities*time_step/space_step).^2; # wave equation coefficient 
ALPHA[nx÷3:nx÷3+10,ny÷3:ny÷3+2] .*= 0.0
ALPHA[2nx÷3:2nx÷3+10,2ny÷3:2ny÷3+2] .*= 0.0
U = zeros(nx, ny); # actual #[exp(- 0.1 * (x-5)^2 - 0.1 * (y-5)^2) for x in xs, y in ys]
V = zeros(nx, ny); # new
# V = [exp(-5*((x-L/3)^2+(y-L/3)^2)) for x in X, y in Y]; # new
V[:,4] .= [exp(-5*x^2) for x in X]
W = zeros(nx, ny); # old

SIM_STEPS = 800;
wave_oh2 = FDTD.oh2(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
wave_oh4 = FDTD.oh4(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
wave_oh6 = FDTD.oh6(U,copy(V),ALPHA, KAPPA, X,Y, SIM_STEPS;_absorbing=true); GC.gc()
end

begin
  GC.gc()
  MAX_VAL = maximum(wave_oh2)
#  sl = SliderGrid(fig[2,1],(label="t", range=2:SIM_STEPS,startvalue=2))
#  hm = lift(sl.sliders[1].value) do _i
  fig = Figure(size=(800,400))
  ax = [Axis(fig[1,1]), Axis(fig[1,2]), Axis(fig[1,3])]
  plot_oh2 = heatmap!(ax[1],X,Y,abs.(wave_oh2[2,:,:]),colormap=:balance,interpolate=true,colorrange=(0,1))
  plot_oh4 = heatmap!(ax[2],X,Y,abs.(wave_oh4[2,:,:]),colormap=:balance,interpolate=true,colorrange=(0,1))
  plot_oh6 = heatmap!(ax[3],X,Y,abs.(wave_oh6[2,:,:]),colormap=:balance,interpolate=true,colorrange=(0,1))
  # vol_data = contour!(ax,X,Y,Z,abs.(volumes[2]), colormap=:turbo, levels=0.1:0.1:0.6)
#  vol_data = contour!(ax,X,Y,Z,volumes[2], transparency=false, colormap=:turbo, levels=0.1:0.1:0.2)
  [ax.aspect=DataAspect() for ax in ax]
  for _i in 3:SIM_STEPS
  #  text!(ax, 0.0, 0.0, text="$_i")
    plot_oh2[3][]=abs.(wave_oh2[_i,:,:])
    plot_oh4[3][]=abs.(wave_oh4[_i,:,:])
    plot_oh6[3][]=abs.(wave_oh6[_i,:,:])
    sleep(0.001)
    println(_i)
    display(fig)
  end
  #volume!(ax, X,Y,Z,volumes[_i])
#  end
  fig
end

# init procedure
begin
L = 5.0
nx,ny,nz   = 51,51,51;
X,Y,Z      = LinRange(-L,L,nx),LinRange(-L,L,ny), LinRange(-L,L,nz);
speed      = 0.2;
space_step = 1.0;# spatial width
time_step  = 1.0; # time step width (remember, dt < dx^2/2)?
velocities = ones(nx,ny,nz) * speed # velocities 
KAPPA      = copy(velocities) * time_step/space_step  # 
ALPHA      = (velocities*time_step/space_step).^2; # wave equation coefficient 
ALPHA[nx÷4:3nx÷4,ny÷4:3ny÷4,nz÷4:3nz÷4] .*= 2.0
U = zeros(nx, ny, nz); # actual #[exp(- 0.1 * (x-5)^2 - 0.1 * (y-5)^2) for x in xs, y in ys]
V = zeros(nx, ny, nz); # new
V = [exp(-2*((x)^2+(y)^2+(z)^2)) for x in X, y in Y, z in Z]; # new
W = zeros(nx, ny, nz); # old

SIM_STEPS = 200;
wave_oh4 = FDTD.oh6(U,V,ALPHA, KAPPA, X,Y,Z, SIM_STEPS;_absorbing=true); GC.gc()
volumes = eachslice(wave_oh4,dims=1)
end

#wave_oh6 = FDTD.oh6(U,V, ALPHA, KAPPA, X,Y,Z, SIM_STEPS); GC.gc()
# wave_oh2 = FDTD.oh2(U,V, ALPHA, KAPPA, X,Y,Z, SIM_STEPS); GC.gc()


begin
  GC.gc()
  MAX_VAL = maximum(wave_oh4)
#  sl = SliderGrid(fig[2,1],(label="t", range=2:SIM_STEPS,startvalue=2))
#  hm = lift(sl.sliders[1].value) do _i
  fig = Figure()
  ax = Axis3(fig[1,1])
  vol_data = volume!(ax,X,Y,Z,abs.(volumes[2]), transparency=true, colormap=:turbo, colorrange=(0,0.75MAX_VAL))
  # vol_data = contour!(ax,X,Y,Z,abs.(volumes[2]), colormap=:turbo, levels=0.1:0.1:0.6)
#  vol_data = contour!(ax,X,Y,Z,volumes[2], transparency=false, colormap=:turbo, levels=0.1:0.1:0.2)
  for _i in 3:SIM_STEPS
  #  text!(ax, 0.0, 0.0, text="$_i")
    vol_data[4][]=volumes[_i]
    sleep(0.001)
    println(_i)
    display(fig)
  end
  #volume!(ax, X,Y,Z,volumes[_i])
#  end
  fig
end
