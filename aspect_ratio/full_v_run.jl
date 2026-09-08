include(relpath("supercritical_elliptical/sup_elliptical.jl","shape_in_wdp"))
include(relpath("subcritical_elliptical/sub_elliptical.jl","shape_in_wdp"))

# different aspect ratios
r = (pi/2)
# we can do an aspect ratio sweep from 0.5 to 1.5
aspects = []
for i in 0.25:0.05:2
    append!(aspects,i)
end
# then find x_len, y_len based off of this 
# x_len/y_len = aspect ratio 
xlens = []
for i in 1:length(0.25:0.05:2)
    append!(xlens, sqrt(r^2/aspects[i]))
end
ylens = []
for i in 1:length(0.25:0.05:2)
    append!(ylens, xlens[i]*aspects[i])
end

### subcritical
sub_vels = 0:0.05:0.95
sub_thrust = zeros(3,length(sub_vels))
sub_eff = zeros(3,length(sub_vels))

Threads.@threads for i in 1:length(sub_vels)
    print("thread ",Threads.threadid(), " is doing v = ",sub_vels[i],"\n")
    A = sub_ellipse(sub_vels[i],2*xlens[6],2*ylens[6],10*xlens[6],10*ylens[6],180,true,false)
    sub_thrust[1,i] = A[11]
    sub_eff[1,i] = A[11]/A[14]
    A = sub_ellipse(sub_vels[i],2*xlens[16],2*ylens[16],10*xlens[16],10*ylens[16],180,true,false)
    sub_thrust[2,i] = A[11]
    sub_eff[2,i] = A[11]/A[14]
    A = sub_ellipse(sub_vels[i],2*xlens[36],2*ylens[36],10*xlens[36],10*ylens[36],180,true,false)
    sub_thrust[3,i] = A[11]
    sub_eff[3,i] = A[11]/A[14]
    print(sub_vels[i],"\n")
end

# supercritical
sup_vels = 1.05:0.025:2
sup_thrust = zeros(3,length(sup_vels))
sup_eff = zeros(3,length(sup_vels))
for i in 1:length(sup_vels)
    # print("thread ",Threads.threadid(), " is doing v = ",sup_vels[i],"\n")
    ### Different sizes of mesh to balance number of discrete points the body
    ### occupies and numerical speed 
    if sup_vels[i] <1.25
        A = sup_opti(2*xlens[6],2*ylens[6],sup_vels[i],350,false)
    elseif 1.25 <= sup_vels[i] <1.5
        A = sup_opti(2*xlens[6],2*ylens[6],sup_vels[i],300,false)
    elseif 1.5 <= sup_vels[i] <1.75
        A = sup_opti(2*xlens[6],2*ylens[6],sup_vels[i],250,false)
    else
        A = sup_opti(2*xlens[6],2*ylens[6],sup_vels[i],180,false)
    end
    sup_thrust[1,i] = A[7]
    sup_eff[1,i] = A[7]/A[13]
    A = sup_opti(2*xlens[16],2*ylens[16],sup_vels[i],250,false)
    sup_thrust[2,i] = A[7]
    sup_eff[2,i] = A[7]/A[13]
    A = sup_opti(2*xlens[36],2*ylens[36],sup_vels[i],180,false,true)
    sup_thrust[3,i] = A[7]
    sup_eff[3,i] = A[7]/A[13]
    print(sup_vels[i],"\n")
end

using JLD2
info = "28/08 full v for divs=250 for r=pi/2, rerunning narrow bodies"
# jldsave("data/full_v_pi_2_div250.jld2";info,sub_vels,sub_thrust,sub_eff,sup_vels,sup_thrust,sup_eff)
# jldsave("full_v_pi_2_narrow.jld2";info,sup_vels,sup_thrust,sup_eff)
### Commented out to prevent overwrite

using Plots
plot(sup_vels,sup_eff[1,:])



