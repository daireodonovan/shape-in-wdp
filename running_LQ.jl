# running LQ checks over the weekend

### Subcritical

include("subcritical_elliptical/sub_elliptical.jl")
using LinearAlgebra
using LaTeXStrings

l_inner_x = pi
l_inner_y = pi

function lq(A,v)
    a = l_inner_x/2
    b = l_inner_y/2
    n = length(A[1][1,:])
    lq_real = zeros(n,n)
    lq_imag = zeros(n,n)
    lq_rcentre = zeros(n,n)
    dx = A[5]
    dy = A[10]
    lhs_x = -5*l_inner_x/2 
    lhs_y = -5*l_inner_y/2
    for xi in 1:n
        for yi in 1:n
            # check if it is in the ellipse
            x = lhs_x + (xi-1)*dx # minus one because we start xi at 1
            y = lhs_y + (yi-1)*dy
            if (x/a)^2 + (y/b)^2 > 1 -2dx
                # do nothing
            elseif (x/a)^2 + (y/b)^2 > 1- 2*dy && (x/a)^2 + (y/b)^2 <= 1 && y<0 && x<0
                re = ((1-v^2)*second_forw_x(A[1],dx,yi,xi)
                        + 2*v*first_forw_x(A[2],dx,yi,xi) 
                        + second_forw_y(A[1],dy,yi,xi)
                        + A[1][yi,xi])
                im = ((1-v^2)*second_forw_x(A[2],dx,yi,xi)
                        - 2*v*first_forw_x(A[1],dx,yi,xi)
                        + second_forw_y(A[2],dy,yi,xi)
                        + A[2][yi,xi])
                lq_real[yi,xi] = abs(re)
                lq_imag[yi,xi] = abs(im)
            elseif (x/a)^2 + (y/b)^2 > 1- 2*dy && (x/a)^2 + (y/b)^2 <= 1 && y>=0 && x>=0
                re = ((1-v^2)*second_back_x(A[1],dx,yi,xi)
                        + 2*v*first_back_x(A[2],dx,yi,xi) 
                        + second_back_y(A[1],dy,yi,xi)
                        + A[1][yi,xi])
                im = ((1-v^2)*second_back_x(A[2],dx,yi,xi)
                        - 2*v*first_back_x(A[1],dx,yi,xi)
                        + second_back_y(A[2],dy,yi,xi)
                        + A[2][yi,xi])
                lq_real[yi,xi] = abs(re)
                lq_imag[yi,xi] = abs(im)
            elseif (x/a)^2 + (y/b)^2 > 1- 2*dy && (x/a)^2 + (y/b)^2 <= 1 && y<0 && x>=0
                re = ((1-v^2)*second_back_x(A[1],dx,yi,xi)
                        + 2*v*first_back_x(A[2],dx,yi,xi) 
                        + second_forw_y(A[1],dy,yi,xi)
                        + A[1][yi,xi])
                im = ((1-v^2)*second_back_x(A[2],dx,yi,xi)
                        - 2*v*first_back_x(A[1],dx,yi,xi)
                        + second_forw_y(A[2],dy,yi,xi)
                        + A[2][yi,xi])
                lq_real[yi,xi] = abs(re)
                lq_imag[yi,xi] = abs(im)
            elseif (x/a)^2 + (y/b)^2 > 1- 2*dy && (x/a)^2 + (y/b)^2 <= 1 && y>=0 && x<0
                re = ((1-v^2)*second_forw_x(A[1],dx,yi,xi)
                        + 2*v*first_forw_x(A[2],dx,yi,xi) 
                        + second_back_y(A[1],dy,yi,xi)
                        + A[1][yi,xi])
                im = ((1-v^2)*second_forw_x(A[2],dx,yi,xi)
                        - 2*v*first_forw_x(A[1],dx,yi,xi)
                        + second_back_y(A[2],dy,yi,xi)
                        + A[2][yi,xi])
                lq_real[yi,xi] = abs(re)
                lq_imag[yi,xi] = abs(im)
            else
                re = ((1-v^2)*second_centred_x(A[1],dx,yi,xi)
                        + 2*v*first_centred_x(A[2],dx,yi,xi) 
                        + second_centred_y(A[1],dy,yi,xi)
                        + A[1][yi,xi])
                im = ((1-v^2)*second_centred_x(A[2],dx,yi,xi)
                        - 2*v*first_centred_x(A[1],dx,yi,xi)
                        + second_centred_y(A[2],dy,yi,xi)
                        + A[2][yi,xi])
                lq_real[yi,xi] = abs(re)
                lq_imag[yi,xi] = abs(im)
                lq_rcentre[yi,xi] = abs(re)
            end
        end
    end
    return norm(lq_real,Inf), norm(lq_imag,Inf),norm(lq_rcentre,Inf),lq_real,lq_imag,lq_rcentre
end

vels = 0:0.025:0.95
max_error_real = zeros(length(vels))
max_error_imag = zeros(length(vels))
max_error_rcentre = zeros(length(vels))

Threads.@threads for i in 1:length(vels)
    print("thread ",Threads.threadid(), " is doing v = ",vels[i],"\n")
    v=vels[i]
    l_inner_x = pi
    l_inner_y = pi
    A = sub_ellipse(v,l_inner_x,l_inner_y,5*l_inner_x,5*l_inner_y,150,true)
    lq_res = lq(A,v)
    
    # store the maximum value
    max_error_real[i] = lq_res[1]
    max_error_imag[i] = lq_res[2]
    max_error_rcentre[i] = lq_res[3]
    print(v,"\t")

end
using JLD2
jldsave("sub_pi_lq_1408.jld2"; vels,max_error_real, max_error_imag,max_error_rcentre)

### Supercritical

include("supercritical_elliptical/sup_elliptical.jl")
theme(:vibrant)
using LinearAlgebra

x_leng=pi
y_leng=pi

function sup_lq(C,v)
    n = div+1
    dx = C[8]
    dy = C[9]
    a = x_leng/2
    b = y_leng/2
    l = x_leng
    dt = dx
    rhs_t = a + dt 
    ymax = round(2*l/sqrt(v^2-1) + l/2,digits=2)
    ymin = round(-2*l/sqrt(v^2-1) - l/2,digits=2)

    # mod_grad_real = zeros(n,length(C[5]))
    # mod_grad_imag = zeros(n,length(C[5]))
    # mod_grad = zeros(n,length(C[5]))
    lq_real = zeros(n,length(C[5]))
    lq_imag = zeros(n,length(C[5]))
    lq_rcentre = zeros(n,length(C[5]))
    for ti in 1:length(C[5])
        for yi in 1:n
            t = -rhs_t + (ti-1)*dt
            y = -ymax + (yi-1)*dy
            if (t/a)^2 + (y/b)^2 > 1 - 2dy
                # nothing
            elseif (t/a)^2 + (y/b)^2 > 1- 3*dy && (t/a)^2 + (y/b)^2 <= 1 && y<0 && t<0
                re = ((1-v^2)*second_forw_x(C[1],dx,yi,ti)
                                - 2*v*first_forw_x(C[2],dx,yi,ti) 
                                + second_forw_y(C[1],dy,yi,ti) 
                                + C[1][yi,ti])
                # re = 1
                im = ((1-v^2)*second_forw_x(C[2],dx,yi,ti) 
                                + 2*v*first_forw_x(C[1],dx,yi,ti)
                                + second_forw_y(C[2],dy,yi,ti) 
                                + C[2][yi,ti])
                lq_real[yi,ti] = abs(re)
                lq_imag[yi,ti] = abs(im)
            elseif (t/a)^2 + (y/b)^2 > 1- 3*dy && (t/a)^2 + (y/b)^2 <= 1 && y>=0 && t>=0
                re = ((1-v^2)*second_back_x(C[1],dx,yi,ti)
                                - 2*v*first_back_x(C[2],dx,yi,ti) 
                                + second_back_y(C[1],dy,yi,ti) 
                                + C[1][yi,ti])
                # re=2
                im = ((1-v^2)*second_back_x(C[2],dx,yi,ti) 
                                + 2*v*first_back_x(C[1],dx,yi,ti)
                                + second_back_y(C[2],dy,yi,ti) 
                                + C[2][yi,ti])
                lq_real[yi,ti] = abs(re)
                lq_imag[yi,ti] = abs(im)
            elseif (t/a)^2 + (y/b)^2 > 1- 3*dy && (t/a)^2 + (y/b)^2 <= 1 && y<0 && t>=0
                re = ((1-v^2)*second_back_x(C[1],dx,yi,ti)
                                - 2*v*first_back_x(C[2],dx,yi,ti) 
                                + second_forw_y(C[1],dy,yi,ti) 
                                + C[1][yi,ti])
                # re = 3
                im = ((1-v^2)*second_back_x(C[2],dx,yi,ti) 
                                + 2*v*first_back_x(C[1],dx,yi,ti)
                                + second_forw_y(C[2],dy,yi,ti) 
                                + C[2][yi,ti])
                lq_real[yi,ti] = abs(re)
                lq_imag[yi,ti] = abs(im)
            elseif (t/a)^2 + (y/b)^2 > 1- 3*dy && (t/a)^2 + (y/b)^2 <= 1 && y>=0 && t<0
                re = ((1-v^2)*second_forw_x(C[1],dx,yi,ti)
                                - 2*v*first_forw_x(C[2],dx,yi,ti) 
                                + second_back_y(C[1],dy,yi,ti) 
                                + C[1][yi,ti])
                im = ((1-v^2)*second_forw_x(C[2],dx,yi,ti) 
                                + 2*v*first_forw_x(C[1],dx,yi,ti)
                                + second_back_y(C[2],dy,yi,ti) 
                                + C[2][yi,ti])
                lq_real[yi,ti] = abs(re)
                lq_imag[yi,ti] = abs(im)
            elseif (t/a)^2 + (y/b)^2 <= 1 - 3*dy
                re = ((1-v^2)*second_centred_x(C[1],dx,yi,ti)
                                - 2*v*first_centred_x(C[2],dx,yi,ti) 
                                + second_centred_y(C[1],dy,yi,ti) 
                                + C[1][yi,ti])
                # re=0
                im = ((1-v^2)*second_centred_x(C[2],dx,yi,ti) 
                                + 2*v*first_centred_x(C[1],dx,yi,ti)
                                + second_centred_y(C[2],dy,yi,ti) 
                                + C[2][yi,ti])
                lq_real[yi,ti] = abs(re)
                lq_imag[yi,ti] = abs(im)
                lq_rcentre[yi,ti] = abs(re)
            end
        end
    end
    return norm(lq_real,Inf), norm(lq_imag,Inf), norm(lq_rcentre,Inf), lq_real, lq_imag, lq_rcentre
end
print("starting supercritical")
div = 100;
re_dif = [];
im_dif = [];
r_cent = [];
x_leng = pi
y_leng = pi
sup_vels = 1.1:0.05:2.5
for i in 1:length(sup_vels)
    # print("thread ",Threads.threadid(), " is doing v = ",sup_vels[i],"\n")
    A = sup_opti(x_leng,y_leng,sup_vels[i],div)
    dif = sup_lq(A,sup_vels[i]);
    append!(re_dif,dif[1])
    append!(im_dif,dif[2])
    append!(r_cent,dif[3])
    print(sup_vels[i],"\n")
end

using JLD2
jldsave("sup_pi_lq_1808.jld2"; sup_vels, re_dif,im_dif,r_cent)
