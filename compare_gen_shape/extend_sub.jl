### double domain size of the subcritical result
using JuMP
import Ipopt
# using SpecialFunctions
include(relpath("functions/integral_functions.jl","shape_in_wdp"))
include(relpath("functions/finite_functions.jl","shape_in_wdp"))
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


function sub_ellipse(v,Qr,Qi,dx,dy,og_len)
    solver=true
    # build optimiser with a div x div grid 
    # solver = true

    divisions = 4*(length(Qr[:,1])-1)
    n = divisions +1;


    # Lorentz factor γ
    gamma = gam(v)

    opti = Model(Ipopt.Optimizer);
    # set_silent(opti);
    ## make our matrices for Q and h
    Q_real = zeros(n,n);
    Q_imag = zeros(n,n);

    @variable(opti, h_real[1:n , 1:n], start=0.01);
    @variable(opti, h_imag[1:n , 1:n], start=0.01);

    # fill in the source matrices
    left_start = Int64(3divisions/8)
    for i in 1:length(Qr[:,1])
        for j in 1:length(Qr[:,1])
            Q_real[i+left_start,j+left_start] = Qr[i,j]
            Q_imag[i+left_start,j+left_start] = Qi[i,j]
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
    l_outer_x=4*og_len
    l_outer_y=4*og_len
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
    @objective(opti, Max, 1);
    optimize!(opti)
    H_real = zeros(n,n);
    H_imag = zeros(n,n);
    for i in 1:n
        for j in 1:n
            H_real[i,j] = value.(h_real[i,j])
            H_imag[i,j] = value.(h_imag[i,j])
        end
    end

    return Q_real, Q_imag, H_real, H_imag
end
