
using JuMP
import Ipopt
# using SpecialFunctions
include(relpath("functions/integral_functions.jl","shape-in-wdp"))
include(relpath("functions/finite_functions.jl","shape-in-wdp"))
using Plots
# using LaTeXStrings

# functions to use

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


function sub_ellipse(v,l_inner_x=3pi/2,l_inner_y=3pi/2,l_outer_x=15pi/2,l_outer_y=15pi/2,divisions=100,silent=false,solver=false)
    # build optimiser with a div x div grid 
    if l_inner_x > l_outer_x || l_inner_y > l_outer_y
        print("inner rectangle has been input as larger than the outer square :(")
        return
    end
    # solver = true

    # divisions = 100;
    n = divisions +1;
    dx = l_outer_x/divisions;
    dy = l_outer_y/divisions;

    # Lorentz factor γ
    gamma = gam(v)

    opti = Model(Ipopt.Optimizer);
    if silent == true
        set_silent(opti);
    end
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

    # for indexing, we have [i,j] which matches to (y,x) because of how a matrix 
    # is row and columns 

    # ==========================
    # Q = 0 ellipse constraint
    # ==========================
    lhs_x = -l_outer_x/2 # the value for the left hand side boundary 
    lhs_y = -l_outer_y/2
    a = l_inner_x/2
    b = l_inner_y/2
    # if we have a set number n = divisions+1
    # also use this loop to find furthest lhs and rhs indices in x and y
    xleft_index = Int64(round(n/2,digits=0))
    xright_index = Int64(round(n/2,digits=0))
    yright_index = Int64(round(n/2,digits=0))
    yleft_index = Int64(round(n/2,digits=0))
    for xi in 1:n
        for yi in 1:n
            # check if it is in the ellipse
            x = lhs_x + (xi-1)*dx # minus one because we start xi at 1
            y = lhs_y + (yi-1)*dy
            if (x/a)^2 + (y/b)^2 > 1
                if solver == false
                    @constraint(opti, Q_real[yi,xi] == 0);
                    @constraint(opti, Q_imag[yi,xi] == 0);
                end
            else
                if solver == true
                    Q_real[yi,xi] = delta_func(dx,lhs_y + (yi-1)*dy,lhs_x + (xi-1)*dx) #0.5*cos.(x/(1+v));
                    # Q_imag[yi,xi] = -0.5*sin.(x/(1+v));
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
                            + second_centred_y(h_real,dy,i,j) 
                            + h_real[i,j] == -Q_real[i,j])
            @constraint(opti, (1-v^2)*second_centred_x(h_imag,dx,i,j) 
                            - 2*v*first_centred_x(h_real,dx,i,j)
                            + second_centred_y(h_imag,dy,i,j) 
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
    lhs_x = -l_outer_x/2
    lhs_y = -l_outer_y/2
    # x- = first_forw_x
    for y_i in 2:n-1
        x_i = 1
        y_val = lhs_y + y_i*dy
        x_val = lhs_x
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_centred_y(h_real,dy,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_centred_y(h_imag,dy,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # x+ = first_back_x
    for y_i in 2:n-1
        x_i = n
        y_val = lhs_y + y_i*dy
        x_val = -lhs_x
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_centred_y(h_real,dy,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_centred_y(h_imag,dy,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # y- = first_forw_y
    for x_i in 2:n-1
        y_i = 1
        y_val = lhs_y
        x_val = lhs_x + x_i*dx 
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dy,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dy,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # y+ = first_back_y
    for x_i in 2:n-1
        y_i = n
        y_val = -lhs_y 
        x_val = lhs_x + x_i*dx 
        #real
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_real,dx,y_i,x_i)
                        + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dy,y_i,x_i)
                        == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
        #imaginary
        @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_centred_x(h_imag,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dy,y_i,x_i)
                    == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    end

    # Corners need to be donw too 
    
    #x-, y-
    x_i = 1 
    y_i = 1
    y_val = lhs_y
    x_val = lhs_x
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dy,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dy,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x-, y+
    x_i = 1 
    y_i = n
    y_val = -lhs_y
    x_val = lhs_x
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dy,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_forw_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dy,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x+, y-
    x_i = n 
    y_i = 1
    y_val = lhs_y
    x_val = -lhs_x
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_forw_y(h_real,dy,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_forw_y(h_imag,dy,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # x+, y+ 
    x_i = n
    y_i = n
    y_val = -lhs_y
    x_val = -lhs_x
    #real
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_real,dx,y_i,x_i)
                    + sin_theta(y_val,x_val,gamma)*first_back_y(h_real,dy,y_i,x_i)
                    == -gamma*h_imag[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));
    #imaginary
    @constraint(opti, (1/gamma)*cos_theta(y_val,x_val,gamma)*first_back_x(h_imag,dx,y_i,x_i)
                + sin_theta(y_val,x_val,gamma)*first_back_y(h_imag,dy,y_i,x_i)
                == gamma*h_real[y_i,x_i]*(1+v*cos_theta(y_val,x_val,gamma)));

    # all boundary conditions are now included 

    # ====================
    # optimisation 
    # ====================
    if solver == true
        @objective(opti, Max, 1);

        l_x = (1/dx)*(l_outer_x/2 - a + dx)
        r_x = (1/dx)*(l_outer_x/2 + a + dx)
        l_y = (1/dy)*(l_outer_y/2 - b + dy)
        r_y = (1/dy)*(l_outer_y/2 + b + dy)

        xright_index = Int64(ceil(r_x)) + 1
        xleft_index = Int64(floor(l_x)) -1
        yleft_index = Int64(floor(l_y)) -1
        yright_index = Int64(ceil(r_y)) + 1

    else
        l_x = (1/dx)*(l_outer_x/2 - a + dx)
        r_x = (1/dx)*(l_outer_x/2 + a + dx)
        l_y = (1/dy)*(l_outer_y/2 - b + dy)
        r_y = (1/dy)*(l_outer_y/2 + b + dy)

        xright_index = Int64(ceil(r_x)) + 1
        xleft_index = Int64(floor(l_x)) -1
        yleft_index = Int64(floor(l_y)) -1
        yright_index = Int64(ceil(r_y)) + 1
        # bounded norm constraint
        @constraint(opti, dx*dy * sum(sum(Q_real[i,j]^2 + Q_imag[i,j]^2 
                for i in yleft_index:yright_index) for j in xleft_index:xright_index) <= 1)#*pi*l_inner_x*l_inner_y/4);

        # maximise the thrust 
        # need to derive the v_thrust 
        @objective(opti, Max,sup_x_thrust_Q(yleft_index, yright_index, xleft_index, xright_index, dx, dy, Q_real, Q_imag,h_real,h_imag))
    end

    optimize!(opti)

    thrust_res=sup_x_thrust_Q(yleft_index, yright_index, xleft_index, xright_index, dx, dy, Q_real, Q_imag,h_real,h_imag)
    pwr = sup_q_pwr(yleft_index, yright_index, xleft_index, xright_index, dx, dy, Q_real,Q_imag,h_real,h_imag,v)
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
    is_optimal = is_solved_and_feasible(opti)
    is_global =  is_solved_and_feasible(opti; allow_local=false)

    return q_real, q_imag, H_real, H_imag, dx, xleft_index, xright_index, yleft_index, yright_index, dy, value.(thrust_res),is_optimal,is_global, value.(pwr)

end
