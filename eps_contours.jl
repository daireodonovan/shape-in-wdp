using JLD2
ENV["MPLBACKEND"] = "Agg"
using Plots
pythonplot()
theme(:vibrant)

# figure 5
data = jldopen("aspect_ratio/data/sub_smaller_2506.jld2")
aspect = data["aspects"];
vels = data["vels"];
thrusts = data["thrusts"];
pwrs = data["pwrs"];
effs = data["effs"];

plts = contour(aspect,vels,thrusts,fill=true,labelfontsize=18,
                legendfontsize=15,tickfontsize=18,margin=0.4Plots.cm,
                xticks=[0.5,1,1.5,2],colorbar_tickfontsize=18)
contour!(aspect,vels,thrusts, c=cgrad([:black,:black]),linewidth=0.5)

savefig(plts, "fig5a.eps")


### figure 8
data = jldopen("general_shape_opti/data/v0_sub_shape.jld2")
D = data["D"];
l = 3pi/2
tic = -5l/2:D[5]:5l/2

r = 3pi/4
circ_x = -3pi/4:pi/128:3pi/4
circ_plus = []
circ_minus = []
for i in 1:length(circ_x)
    append!(circ_plus, sqrt(r^2 - circ_x[i]^2))
    append!(circ_minus, -sqrt(r^2 - circ_x[i]^2))
end

theme(:default)
plh = heatmap(tic,tic,D[3],tickfontsize=18,
        clims=(-0.31,0.3),size=(400,400),colorbar=false)
plot!(plh,circ_x,circ_plus,label=false,color=:lightgreen,linestyle=:dash,
        linewidth=3)
plot!(plh,circ_x,circ_minus,label=false,color=:lightgreen,linestyle=:dash,  
        linewidth=3)

savefig(plh,"fig8a_h.eps")

hcbar = heatmap(rand(2,2), framestyle=:none, cbar=true, lims=(-1,0),
            size=(500,400),tickfontsize=16,right_margin=1.25Plots.cm,
            clims=(-0.31,0.31),colorbar_tickfontsize=18)
savefig(hcbar,"fig8_hbar.eps")

qcbar = heatmap(rand(2,2), framestyle=:none, cbar=true, lims=(-1,0),
            size=(500,200),tickfontsize=16,right_margin=1.25Plots.cm,
            clims=(-0.6,0.6),c=:seismic,colorbar_ticks=[0],colorbar_tickfontsize=0)
savefig(qcbar,"fig8_qbar.eps")

plt_qb = heatmap(tic,tic,D[1],clims=(-0.6,0.6),c=:seismic,
        xlims=[-3pi/3,3pi/3],ylims=[-3pi/3,3pi/3],
        tickfontsize=16,
        size=(400,400),xticks=false,yticks=false,colorbar=false)
plot!(plt_qb, circ_x,circ_plus,label=false,color="#00ff22",linestyle=:dot,
        linewidth=3)
plot!(plt_qb, circ_x,circ_minus,label=false,color="#00ff22",linestyle=:dot,  
        linewidth=3)
savefig(plt_qb,"fig8a_q.eps")

data5 = jldopen("general_shape_opti/data/v05_sub_shape.jld2")
A = data5["A"];
B = data5["B"];
tic2 = -5l/2:A[5]:5l/2
plt_q2b = heatmap(tic2,tic2,A[1],clims=(-0.6,0.6),c=:seismic,
        xlims=[-3pi/3,3pi/3],ylims=[-3pi/3,3pi/3],
        tickfontsize=16,
        size=(400,400),xticks=false,yticks=false,colorbar=false)
plot!(plt_q2b, circ_x,circ_plus,label=false,color="#00ff22",linestyle=:dot,
        linewidth=3)
plot!(plt_q2b, circ_x,circ_minus,label=false,color="#00ff22",linestyle=:dot,  
        linewidth=3)
savefig(plt_q2b, "fig8b_q.eps")

plh = heatmap(tic2,tic2,A[3],tickfontsize=15,
        clims=(-0.31,0.3),size=(400,400),colorbar=false)
plot!(plh,circ_x,circ_plus,label=false,color=:lightgreen,linestyle=:dash,
        linewidth=3)
plot!(plh,circ_x,circ_minus,label=false,color=:lightgreen,linestyle=:dash,  
        linewidth=3)
savefig(plh,"fig8b_h.eps")

### fig 7 - make more data with this 
using LaTeXStrings
data = jldopen("aspect_ratio/data/full_v_sweep.jld2")
sub_vels = data["sub_vels"];
sup_vels = data["sup_vels"];
sub_thrust = data["sub_thrust"];
sup_thrust = data["sup_thrust"];
y = ones(5)
x = 0.9:0.05:1.1
length(x)
plt = plot(x,y,label=false,fillrange=0,color="grey", fillstyle = :/)
plot!(plt,sub_vels,sub_thrust[1,:],
    linewidth=3,
    labelfontsize=16,tickfontsize=14,legendfontsize=16,gridalpha=0.5
    ,ylims=[0,0.5],color="#0072B2",label=L"\mathrm{aspect}=0.5",
    size=(800,300))
plot!(plt,sup_vels,new_sup[1,:],linewidth=3,color="#0072B2",label=false)
plot!(plt,sub_vels,sub_thrust[2,:],
    linewidth=3,
    color="#DC267F",label=L"\mathrm{aspect}=1")
plot!(plt,sup_vels,new_sup[2,:],linewidth=3,color="#DC267F",label=false)
plot!(plt,sub_vels,sub_thrust[3,:],
    linewidth=3,
    color="#E69F00",label=L"\mathrm{aspect}=2")
plot!(plt,sup_vels,new_sup[3,:],linewidth=3,color="#E69F00",label=false)
