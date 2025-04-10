using Images
using GLMakie; GLMakie.activate!(inline=false; float=true)

coarse = Float64.(Gray.(load("./images/coarse_crop.jpg")))
fine   = Float64.(Gray.(load("./images/fine_crop.jpg")))

edges_fine = reverse(size(fine))
edges_coarse = reverse(size(coarse))


let
_cmap = :turku
fig = Figure(size = (3,3) .* 200)
ax1 = Makie.Axis(fig[1, 1])
ax2 = Makie.Axis(fig[2, 1])
ax_inset = Makie.Axis(fig[1,1],
                      width=Relative(0.2),
                      height=Relative(0.2),
                      halign=0.88,
                      valign=0.2,)
heatmap!(ax1, rotr90(fine),        colormap=_cmap,)# colorrange=(0,1))
hm=heatmap!(ax2, rotr90(coarse),   colormap=_cmap,)# colorrange=(0,1))
heatmap!(ax_inset, rotr90(coarse), colormap=_cmap,)# colorrange=(0,1))
hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax_inset)
ax1.aspect=DataAspect()
ax2.aspect=DataAspect()
ax_inset.aspect=DataAspect()
# Label(fig[1,0], "fine", rotation=pi/2, tellheight=false)

w = 290
w2 = 240

x1 = 3700
x2 = x1 + w2
y1 = 585
y2 = y1 + w2

# poly!(ax1, Point2f[(0,edges_fine[2]), (0, edges_fine[2]-w), (w, edges_fine[2]-w), (w,edges_fine[2])], color=:black)
# poly!(ax2, Point2f[(0,edges_coarse[2]), (0, edges_coarse[2]-w), (w, edges_coarse[2]-w), (w,edges_coarse[2])], color=:black)
# poly!(ax1, Point2f[(x1,y1), (x1, y2), (x2, y2), (x2,y1)], color=:black)
# # scatter!(ax1, Point2f(115,2100), color="black", alpha=0.5, marker=:rect,markersize=58)
# # scatter!(ax2, Point2f(115,2100), color="black", alpha=0.5, marker=:rect,markersize=58)
# # scatter!(ax1, )
# text!(ax1, 30, 2050, text="(a)", fontsize = 20, color="white")
# text!(ax2, 30, 2100, text="(b)", fontsize = 20, color="white")
# text!(ax1, x1+20, y1+30, text="(b)", fontsize = 20, color="white")
# cb=Colorbar(fig[1:2,2], hm)
# cb.ticks=[minimum(coarse), maximum(coarse)]

ax_inset.bottomspinecolor=(:black, 0.8)
ax_inset.topspinecolor=(:black, 0.8)
ax_inset.leftspinecolor=(:black, 0.8)
ax_inset.rightspinecolor=(:black, 0.8)
save("/Users/ar/Desktop/coarse_graining_$_cmap.jpg", fig, px_per_unit = 4)

# resize_to_layout!(fig.layout)

fig
end