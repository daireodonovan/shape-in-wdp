# High density grid for paper figures and results
using JLD2
l = pi
shape_l = sqrt(3)*pi # because of the shape being one third of the area

### subcritical 
print("----------------\nstarting\n--------------\n")
include(relpath("subcritical_elliptical/sub_elliptical.jl","shape-in-wdp"))

A = sub_ellipse(0,l,l,5l,5l,240)

B = sub_ellipse(0.5,l,l,5l,5l,240)

info = "data gathering for v=0,v=0.5 for radius pi/2"
# jldsave("sub_pi_0708.jld2";info,A,B) # commented out to prevent overwrite
print("sub done\n")

#### subcritical optimal shaping 

include(relpath("general_shape_opti/sub_gen_shape.jl","shape-in-wdp"))


C = sub_ellipse(0,1,shape_l,shape_l,5*shape_l,5*shape_l,240)

D = sub_ellipse(0.5,1,shape_l,shape_l,5*shape_l,5*shape_l,240)

info = "data gathering for optimal shape v=0,v=0.5 for radius pi/2"
# jldsave("sub_shape_pi_0708.jld2";info,C,D)
print("sub shape done\n")

#### supercritical 
include(relpath("supercritical_elliptical/sup_elliptical.jl","shape-in-wdp"))


E = sup_opti(l, l, 1.25, 250,true)

F = sup_opti(l, l, 1.75, 250,true)

info = "data gathering for v=1,25,v=1.75 for radius pi/2"
# jldsave("sup_pi_0608.jld2";info,E,F)
print("sup done\n")

#### supercritical optimal shaping 
include(relpath("general_shape_opti/sup_gen_shape.jl","shape-in-wdp"))


G = input_opti(shape_l, shape_l, 1.25, 250)

H = input_opti(shape_l, shape_l, 1.75, 250)

info = "data gathering for optimal shape v=1,25,v=1.75 for radius pi/2"
# jldsave("sup_shape_pi_0608.jld2";info,G,H)
print("sup shape done\n")