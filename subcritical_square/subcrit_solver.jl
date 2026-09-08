# Storing the numerical solver somewhere so that I do not have to copy and paste
# it into each file

using JuMP
import Ipopt
using SpecialFunctions
include(relpath("functions/integral_functions.jl","shape-in-wdp"))
include(relpath("functions/finite_functions.jl","shape-in-wdp"))
using Plots
using LaTeXStrings

# main functions in use

"""
    Defines the 2d regularised Dirac delta distribution.
"""
function delta_func(a,y,x)
    return (1/(2*pi*(a^2)))*exp(-(x^2+y^2)/(2*a^2))
end

"""
    Calculates the Lorentz factor γ
"""
function gam(vel)
    return (1-vel^2)^(-1/2)
end


"""
    cos(θ) effective angle
"""
function cos_theta(y,x,gamma)
    return (x*gamma)/(sqrt(x^2*gamma^2 + y^2))
end

"""
    sin(θ) effective angle
"""
function sin_theta(y,x,gamma)
    return y/(sqrt(x^2*gamma^2 + y^2))
end


function vel_opti(v,l_inner=3pi/2,l_outer=15pi/2,divisions=100,solver=false)
    # build optimiser with a 100x100 grid 
    if l_inner > l_outer
        print("inner square has been input as larger than the outer square :(")
        return
    end
    

    # divisions = 100;
    n = divisions +1;
    dx = l_outer/divisions;

    # Lorentz factor γ
    gamma = gam(v)

    opti = Model(Ipopt.Optimizer);
    # set_silent(opti);

    ## make our matrices for Q and h
    if solver == false
        @variable(opti, Q_real[1:n , 1:n], start=0.1);
        @variable(opti, Q_imag[1:n , 1:n], start=0.1);
    else
        Q_real = zeros(n,n);
        Q_imag = zeros(n,n);
    end

    @variable(opti, h_real[1:n , 1:n], start=0.01);
    @variable(opti, h_imag[1:n , 1:n], start=0.01);
    # find square
    if l_inner == l_outer
        left = "N/A"
        right="N/A"
        left_index = 1
        right_index = n
    else
        left = nothing
        right = nothing
        left_index=0
        right_index=0
        for j in 1:n
            if -l_inner/2 <= -l_outer/2 + j*dx && left == nothing
                if abs(abs(-l_outer/2 + (j-1)*dx)-3pi/4) - abs(abs(-l_outer/2 + j*dx)-3pi/4) < 0
                    left = abs(-l_outer/2 + (j-1)*dx) - 3pi/4
                    left_index = j-1
                else
                    left = abs(-l_outer/2 + j*dx)-3pi/4
                    left_index = j
                end
            end

            if -l_outer/2 + j*dx > l_inner/2 && right == nothing
                if abs(abs(-l_outer/2 + (j-1)*dx)-3pi/4) - abs(abs(-l_outer/2 + j*dx)-3pi/4) < 0
                    right = abs(-l_outer/2 + (j-1)*dx) - 3pi/4
                    right_index = j-1
                else
                    right = abs(-l_outer/2 + j*dx)-3pi/4
                    right_index = j
                end
            end
        end
        raft_indices = 0
        for i in 1:n
            for j in 1:n
                if left_index <= i <= right_index && left_index <= j <= right_index
                    raft_indices += 1
                    if solver == true
                        lhs = -l_outer/2
                        Q_real[i,j] = delta_func(dx,lhs + (i-1)*dx,lhs + (j-1)*dx);# 0.1*cos.((lhs + j*dx)/(1+v));
                        # Q_imag[i,j] = -0.1*sin.((lhs + j*dx)/(1+v));
                    end
                else
                    if solver == false
                        @constraint(opti, Q_real[i,j] == 0);
                        @constraint(opti, Q_imag[i,j] == 0);
                    end
                end
            end
        end
    end

    # ====================
    # define wave equation
    # ====================
    for i in 2:n-1
        for j in 2:n-1
            @constraint(opti, (1-v^2)*second_centred_x(h_real,dx,i,j)
                            + 2*v*first_centred_x(h_imag,dx,i,j) 
                            + second_centred_y(h_real,dx,i,j) 
                            + h_real[i,j] == -Q_real[i,j])
            @constraint(opti, (1-v^2)*second_centred_x(h_imag,dx,i,j) 
                            - 2*v*first_centred_x(h_real,dx,i,j)
                            + second_centred_y(h_imag,dx,i,j) 
                            + h_imag[i,j] == -Q_imag[i,j])
        end
    end
    
    # ==========================
    # define boundary conditions
    # ==========================
    # because we have at the edges, we will need to divide the prescription into
    # parts and use forward and backward difference where necessary 

    # The set-up is as follows
    #           y+
        ##################
        #                #
        #                #
    #x- #                # x+
        #                #
        #                #
        #                #
        ##################
    #            y-

    # The boundary condition needs to be divided into real and imaginary parts
    lhs = -l_outer/2
    # x- = first_forw_x
    for y_i in 2:n-1
        x_i = 1
        y_val = lhs + y_i*dx
        x_val = lhs
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_centred_y(h_real,dx,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_centred_y(h_imag,dx,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # x+ = first_back_x
    for y_i in 2:n-1
        x_i = n
        y_val = lhs + y_i*dx
        x_val = -lhs
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_centred_y(h_real,dx,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_centred_y(h_imag,dx,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # y- = first_forw_y
    for x_i in 2:n-1
        y_i = 1
        y_val = lhs 
        x_val = lhs + x_i*dx 
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dx,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dx,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # y+ = first_back_y
    for x_i in 2:n-1
        y_i = n
        y_val = -lhs 
        x_val = lhs + x_i*dx 
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dx,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dx,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # Corners need to be donw too 
    
    #x-, y-
    x_i = 1 
    y_i = 1
    y_val = lhs 
    x_val = lhs 
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dx,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dx,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x-, y+
    x_i = 1 
    y_i = n
    y_val = -lhs 
    x_val = lhs 
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dx,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dx,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x+, y-
    x_i = n 
    y_i = 1
    y_val = lhs 
    x_val = -lhs 
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dx,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dx,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x+, y+ 
    x_i = n
    y_i = n
    y_val = -lhs 
    x_val = -lhs 
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dx,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dx,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));


    # all boundary conditions are now included 

    # ====================
    # optimisation 
    # ====================
    if solver == true
        @objective(opti, Max, 1);

        optimize!(opti)

        # extract the results
        q_real = zeros(n,n);
        q_imag = zeros(n,n);
        H_real = zeros(n,n);
        H_imag = zeros(n,n);
        for i in 1:n
            for j in 1:n
                q_real[i,j] = value.(Q_real[i,j])
                q_imag[i,j] = value.(Q_imag[i,j])
                H_real[i,j] = value.(h_real[i,j])
                H_imag[i,j] = value.(h_imag[i,j])
            end
        end
        thrust_res = x_thrust_Q(left_index, right_index, left_index, right_index, dx, value.(Q_real), value.(Q_imag),value.(h_real),value.(h_imag))

        return q_real, q_imag, H_real, H_imag, dx, thrust_res
    else
        # bounded norm constraint
        @constraint(opti, dx^2 * sum(sum(Q_real[i,j]^2 + Q_imag[i,j]^2 
                for i in left_index-1:right_index+1) for j in left_index-1:right_index+1) <= 1);
        # maximise the thrust 
        # need to derive the v_thrust 
        @objective(opti, Max,x_thrust_Q(left_index, right_index, left_index, right_index, dx, Q_real, Q_imag,h_real,h_imag))

        optimize!(opti)

        # extract the results
        q_real = zeros(n,n);
        q_imag = zeros(n,n);
        H_real = zeros(n,n);
        H_imag = zeros(n,n);
        for i in 1:n
            for j in 1:n
                q_real[i,j] = value.(Q_real[i,j])
                q_imag[i,j] = value.(Q_imag[i,j])
                H_real[i,j] = value.(h_real[i,j])
                H_imag[i,j] = value.(h_imag[i,j])
            end
        end
        thrust_res = x_thrust_Q(left_index, right_index, left_index, right_index, dx, value.(Q_real), value.(Q_imag),value.(h_real),value.(h_imag))


        return q_real, q_imag, H_real, H_imag, dx, left_index, right_index, thrust_res


    end


end

function LQ_check(A,n,v)
    dx = A[5]
    # run the operator L check
    check_real = zeros(n,n);
    check_imag = zeros(n,n);
    for i in 2:n-1
    	for j in 2:n-1
       		check_real[i,j] =(1-v^2)*second_centred_x(A[1],dx,i,j)
                            + 2*v*first_centred_x(A[2],dx,i,j) 
                            + second_centred_y(A[1],dx,i,j) 
                            + A[1][i,j]
        	check_imag[i,j] =(1-v^2)*second_centred_x(A[2],dx,i,j) 
                            - 2*v*first_centred_x(A[1],dx,i,j)
                            + second_centred_y(A[2],dx,i,j) 
                            + A[2][i,j] 
	    end
    end	
    return check_real, check_imag
end

# going to try and do the other boundary conditions and for that need a proper
# angle check
function angle_check(y,x,divs)
    if x == divs/2 && y < divs/2
        return -pi/2
    elseif x == divs/2 && y > divs/2
        return pi/2
    # elseif x == 1 && y== divs/2
    #     return pi # works even without this one
    else
        return atan.((y-divs/2),(x-divs/2)) 
        # modified to work with flipped matrix
    end
end


