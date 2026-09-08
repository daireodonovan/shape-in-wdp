# want to do an aspect ratio sweep 
include(relpath("supercritical_elliptical/sup_elliptical.jl","shape_in_wdp"))
include(relpath("subcritical_elliptical/sub_elliptical.jl","shape_in_wdp"))

### Code to do a sweep 
# want to maintain constant area, so do xlens and ylens so that πab=const

# want a linear relationship in aspect ratio for the x axis

# given radius r 
r = (3pi/4)
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


vstep = 0.05
startv = 1.1
endv = 2.5
# lstep = pi/12

effs = zeros(length(startv:vstep:endv),length(xlens))
thrusts = zeros(length(startv:vstep:endv),length(xlens))
pwrs = zeros(length(startv:vstep:endv),length(xlens))
vels = []
for i in startv:vstep:endv
    append!(vels,i)
end

# uses multiple threads to spead up sweep
Threads.@threads for vi in 1:length(startv:vstep:endv)
    print("thread ",Threads.threadid(), " is doing v = ",vels[vi],"\n")
    for l in 1:length(aspects)
        # multiply the length by 2 so that since currently we are working with half major and minor axes

        #### Subcritical 
        A = sub_ellipse(vels[vi],2*xlens[l],2*ylens[l],xlens[l]*10,ylens[l]*10,100,true) 
        thrust = sup_x_thrust_Q(A[8], A[9], A[6], A[7], A[5], A[10], A[1], A[2],A[3],A[4])
        power = sup_q_pwr(A[8], A[9], A[6], A[7], A[5], A[10],A[1], A[2],A[3],A[4],vels[vi])

        #### Supercritical 
        # A = sup_opti(2*xlens[l],2*ylens[l],vels[vi],150,false,true)
        # power = A[13]      
        # thrust = A[7]
        eff = thrust/power
        effs[vi,l] = eff
        thrusts[vi,l] = thrust
        pwrs[vi,l] = power
        # print("\t")
    end
    print(vels[vi]," ,thread= ", Threads.threadid(),"\t")
end
