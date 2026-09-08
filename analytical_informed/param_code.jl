# code to do analytically informed optimisation

# imports
using JuMP
import Ipopt
using SpecialFunctions
include(relpath("functions/integral_functions.jl","shape_in_wdp"))
include(relpath("functions/finite_functions.jl","shape_in_wdp"))
using Plots
using LaTeXStrings


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


function optimise_param(l_inner, l_outer, divisions,print_values=false)
    # cancels the use of the if the inner square is larger than the outer square
    if l_inner > l_outer
        print("inner square has been input as larger than the outer square :(")
        return
    end
    
    # If all is above board, we will run the function like we normally do
    n = divisions+1
    dx = l_outer/divisions

    opti = Model(Ipopt.Optimizer);
    # set_silent(opti)

    #define the variables as matrices as opposed to vectors.
    # putting start values close to one set of optimal values to see if it 
    # speeds things up.
    @variable(opti, Q_real[1:n , 1:n], start=0.1);
    @variable(opti, Q_imag[1:n , 1:n], start=0.1);
    @variable(opti, A_real,start=0.1);
    @variable(opti, A_imag,start=0.1);
    @variable(opti, B_real,start=0.1);
    @variable(opti, B_imag,start=0.1);
    @variable(opti, C_real,start=0.1);
    @variable(opti, C_imag,start=0.1);
    @variable(opti, D_real,start=0.1);
    @variable(opti, D_imag,start=0.1);
    @variable(opti, k_num<=1,start=0.1);    # do two cases where it is < or >1 


    @variable(opti, h_real[1:n , 1:n], start=0.01);
    @variable(opti, h_imag[1:n , 1:n], start=0.01);

    # for indexing, we have [i,j] which matches to (y,x) because of how a matrix 
    # is row and columns 

    # adding in some Q = 0 constraints in the event we're working in the 
    # extended domain

    # need to specify the values where the body is non-zero with l_inner

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
    end
    for i in 1:n
        for j in 1:n
            if left_index <= i <= right_index && left_index <= j <= right_index
                y_val = -l_outer/2 + dx*i
                x_val = -l_outer/2 + dx*j
                @constraint(opti, Q_real[i,j] == 
                            A_real*cos(k_num*x_val + y_val*sqrt(1-k_num^2))
                            - A_imag*sin(k_num*x_val + y_val*sqrt(1-k_num^2))
                            + B_real*cos(-k_num*x_val + y_val*sqrt(1-k_num^2))
                            - B_imag*sin(-k_num*x_val + y_val*sqrt(1-k_num^2))
                            + C_real*cos(k_num*x_val - y_val*sqrt(1-k_num^2))
                            - C_imag*sin(k_num*x_val - y_val*sqrt(1-k_num^2))
                            + D_real*cos(-k_num*x_val - y_val*sqrt(1-k_num^2))
                            - D_imag*sin(-k_num*x_val - y_val*sqrt(1-k_num^2)) )
                @constraint(opti, Q_imag[i,j] == 
                            A_imag*cos(k_num*x_val + y_val*sqrt(1-k_num^2))
                            + A_real*sin(k_num*x_val + y_val*sqrt(1-k_num^2))
                            + B_real*sin(-k_num*x_val + y_val*sqrt(1-k_num^2))
                            + B_imag*cos(-k_num*x_val + y_val*sqrt(1-k_num^2))
                            + C_real*sin(k_num*x_val - y_val*sqrt(1-k_num^2))
                            + C_imag*cos(k_num*x_val - y_val*sqrt(1-k_num^2))
                            + D_real*sin(-k_num*x_val - y_val*sqrt(1-k_num^2))
                            + D_imag*cos(-k_num*x_val - y_val*sqrt(1-k_num^2)) )
                
            else
                @constraint(opti, Q_real[i,j] == 0);
                @constraint(opti, Q_imag[i,j] == 0);
            end
        end
    end

    for i in 2:n-1
        for j in 2:n-1
            @constraint(opti, second_centred_x(h_real,dx,i,j)
                            + second_centred_y(h_real,dx,i,j) 
                            + h_real[i,j] == -Q_real[i,j])
            @constraint(opti, second_centred_x(h_imag,dx,i,j) 
                            + second_centred_y(h_imag,dx,i,j) 
                            + h_imag[i,j] == -Q_imag[i,j])
        end
    end

    # Had to deconstruct the constraints into loops to use the angle finder
    for k in 2:n-1
        # first check angle
        theta = angle_check(k,1,divisions)
        # then add constraint
        #x- real 
        @constraint(opti, cos.(theta)*first_forw_x(h_real,dx,k,1) + sin.(theta)*first_centred_y(h_real,dx,k,1) == -h_imag[k,1]);
        # This  behaves with JuMP
        #x- imag
        @constraint(opti, cos.(theta)*first_forw_x(h_imag,dx,k,1) + sin.(theta)*first_centred_y(h_imag,dx,k,1) == h_real[k,1]); 
    end

    for k in 2:n-1
        # first check angle
        theta = angle_check(k,n,divisions)
        # then add constraint
        #x+ real 
        @constraint(opti, cos.(theta)*first_back_x(h_real,dx,k,n) + sin.(theta)*first_centred_y(h_real,dx,k,n) == -h_imag[k,n]);
        #x+ imag;
        @constraint(opti, cos.(theta)*first_back_x(h_imag,dx,k,n) + sin.(theta)*first_centred_y(h_imag,dx,k,n) == h_real[k,n]); 
    end

    for k in 2:n-1
        # first check angle
        theta = angle_check(n,k,divisions)
        # then add constraint
        #y- real 
        @constraint(opti, cos.(theta)*first_centred_x(h_real,dx,n,k) + sin.(theta)*first_back_y(h_real,dx,n,k) == -h_imag[n,k]);
        #y- imag
        @constraint(opti, cos.(theta)*first_centred_x(h_imag,dx,n,k) + sin.(theta)*first_back_y(h_imag,dx,n,k) == h_real[n,k]); 
        # swapped sign to see if it fixes
    end

    for k in 2:n-1
        # first check angle
        theta = angle_check(1,k,divisions)
        # then add constraint
        #y+ real 
        @constraint(opti, cos.(theta)*first_centred_x(h_real,dx,1,k) + sin.(theta)*first_forw_y(h_real,dx,1,k) == -h_imag[1,k]);
        #y+ imag
        @constraint(opti, cos.(theta)*first_centred_x(h_imag,dx,1,k) + sin.(theta)*first_forw_y(h_imag,dx,1,k) == h_real[1,k]);
        # swapped sign to see if it fixes
    end


    #############
    ## Corners ##
    #############

    #= 
        x+,y+ = pi/4
        x+,y- = -pi/4
        x-,y+ = 3pi/4
        x-,y- = -3pi/4
    
        I used the angle_check() function instead of precribing them 
    =#
    theta = angle_check(1,n,divisions)
    # x+,y+
    @constraint(opti, cos.(theta)*first_back_x(h_real,dx,1,n) + sin.(theta)*first_forw_y(h_real,dx,1,n) == -h_imag[1,n]);
    @constraint(opti, cos.(theta)*first_back_x(h_imag,dx,1,n) + sin.(theta)*first_forw_y(h_imag,dx,1,n) == h_real[1,n]); 

    # x+,y-
    theta = angle_check(n,n,divisions)
    @constraint(opti, cos.(theta)*first_back_x(h_real,dx,n,n) + sin.(theta)*first_back_y(h_real,dx,n,n) == -h_imag[n,n]);
    @constraint(opti, cos.(theta)*first_back_x(h_imag,dx,n,n) + sin.(theta)*first_back_y(h_imag,dx,n,n) == h_real[n,n]); 

    # x-,y+
    theta = angle_check(1,1,divisions)
    @constraint(opti, cos.(theta)*first_forw_x(h_real,dx,1,1) + sin.(theta)*first_forw_y(h_real,dx,1,1) == -h_imag[1,1]);
    @constraint(opti, cos.(theta)*first_forw_x(h_imag,dx,1,1) + sin.(theta)*first_forw_y(h_imag,dx,1,1) == h_real[1,1]); 

    # x-,y-
    theta = angle_check(n,1,divisions)
    @constraint(opti, cos.(theta)*first_forw_x(h_real,dx,n,1) + sin.(theta)*first_back_y(h_real,dx,n,1) == -h_imag[n,1]);
    @constraint(opti, cos.(theta)*first_forw_x(h_imag,dx,n,1) + sin.(theta)*first_back_y(h_imag,dx,n,1) == h_real[n,1]); 
            # also flipped signs of this (needs checking)

    # still has the correct number of constraints


    ##################
    ## bounded norm ##
    ##################

    @constraint(opti, dx^2 * sum(sum(Q_real[i,j]^2 
        + Q_imag[i,j]^2 for i in left_index:right_index) for j in left_index:right_index) <= 1);


    # thrust with just the fore-aft
    @objective(opti, Max,x_thrust_Q(left_index, right_index, left_index, right_index, dx, Q_real, Q_imag,h_real,h_imag))

    optimize!(opti)

    if print_values == true

    nrm = dx^2 *sum(sum(value.(value.(value.(value.(value.(Q_real)))))[i,j]^2 + value.(Q_imag)[i,j]^2 for i in left_index:right_index) for j in left_index:right_index)
    end
    thrust_res = x_thrust_Q(left_index, right_index, left_index, right_index, dx, value.(Q_real), value.(Q_imag),value.(h_real),value.(h_imag))
    pwr_res = 2*dx^2 * sum(sum(-value.(Q_real)[i,j]*value.(h_imag)[i,j] +
                        value.(Q_imag)[i,j]*value.(h_real)[i,j]  for i in left_index:right_index) for j in left_index:right_index)


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

    # make an array of the output values
    values_to_output = [value.(A_real),value.(A_imag),value.(B_real),value.(B_imag),value.(C_real),value.(C_imag),
                        value.(D_real),value.(D_imag),value.(k_num)]

    return q_real, q_imag, H_real, H_imag, dx, left_index, right_index, values_to_output, thrust_res, pwr_res
end
