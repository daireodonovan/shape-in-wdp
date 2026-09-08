# higher order derivative optimisation code 

include("finite_functions.jl")
include("integral_functions.jl")

using JuMP
import Ipopt

"Function to recover the angle for radial boundary conditions"
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


"""
    Radially symmetric boundary conditions being used with higher order finite
    difference approximations to hopefully show that it gets closer to the 
    semi-analytical optimisation.
"""
function optimise_radial_bc_ho(l_inner, l_outer, divisions,print_values=false,solver=false,nf=false)
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
    set_time_limit_sec(opti,360) # times out after 6 minutes
    #define the variables as matrices as opposed to vectors
    if solver == false
        @variable(opti, Q_real[1:n , 1:n], start=0.1);
        @variable(opti, Q_imag[1:n , 1:n], start=0.1);
    else
        Q_real = zeros(n,n);
        Q_imag = zeros(n,n);
    end

    # Q_real = zeros(n,n)
    # Q_imag = zeros(n,n)

    @variable(opti, h_real[1:n , 1:n], start=0.01);
    @variable(opti, h_imag[1:n , 1:n], start=0.01);

    # for indexing, we have [i,j] which matches to (y,x) because of how a matrix 
    # is row and columns 

    # adding in some Q = 0 constraints in the event we're working in the 
    # extended domain

    # need to specify the values where the body is non-zero with l_inner

    ########################
    ## I overdid this lol ##
    ########################
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

    # counter = 0 
    # A = zeros(n,n)
    for i in 1:n
        for j in 1:n
            
            if left_index <= i <= right_index && left_index <= j <= right_index
                # counter += 1
                # A[i,j] = 1
                # prescribing Q
                if solver == true
                    Q_real[i,j] = 1;
                    # Q_imag[i,j] = 1;
                    # Q_real[i,j] = cos.((-l_outer/2+j*dx))
                    # Q_imag[i,j] = sin.((-l_outer/2+j*dx))
                    # Q_real[i,j] = sin.(sqrt((-l_outer/2+j*dx)^2 + (l_outer/2 - i*dx)^2))
                end
            else
                if solver == false
                    @constraint(opti, Q_real[i,j] == 0);
                    @constraint(opti, Q_imag[i,j] == 0);
                end
            end
        end
    end
    # print("\n","counter = ", counter)
    
    # Because we have higher order finite differences, we will have to do the 
    # edges with forward and backward difs
    for i in 5:n-4
        for j in 5:n-4
            @constraint(opti, second_centred_x_ho4(h_real,dx,i,j) +second_centred_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_centred_x_ho4(h_imag,dx,i,j) +second_centred_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end

    # helmholtz at the edges
    # bottom 
    for i in 2:4
        for j in 5:n-4
            @constraint(opti, second_centred_x_ho4(h_real,dx,i,j) +second_forw_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_centred_x_ho4(h_imag,dx,i,j) +second_forw_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end
    # top 
    for i in n-3:n-1
        for j in 5:n-4
            @constraint(opti, second_centred_x_ho4(h_real,dx,i,j) +second_back_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_centred_x_ho4(h_imag,dx,i,j) +second_back_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end
    # left
    for i in 5:n-4
        for j in 2:4
            @constraint(opti, second_forw_x_ho4(h_real,dx,i,j) +second_centred_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_forw_x_ho4(h_imag,dx,i,j) +second_centred_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end  
    # right
    for i in 5:n-4
        for j in n-3:n-1
            @constraint(opti, second_back_x_ho4(h_real,dx,i,j) +second_centred_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_back_x_ho4(h_imag,dx,i,j) +second_centred_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end  
    # top right 
    for i in n-3:n-1
        for j in n-3:n-1
            @constraint(opti, second_back_x_ho4(h_real,dx,i,j) +second_back_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_back_x_ho4(h_imag,dx,i,j) +second_back_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end      
    # top left
    for i in n-3:n-1
        for j in 2:4
            @constraint(opti, second_forw_x_ho4(h_real,dx,i,j) +second_back_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_forw_x_ho4(h_imag,dx,i,j) +second_back_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end            
    # bottom left
    for i in 2:4
        for j in 2:4
            @constraint(opti, second_forw_x_ho4(h_real,dx,i,j) +second_forw_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_forw_x_ho4(h_imag,dx,i,j) +second_forw_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end  
    # bottom right
    for i in 2:4
        for j in n-3:n-1
            @constraint(opti, second_back_x_ho4(h_real,dx,i,j) +second_forw_y_ho4(h_real,dx,i,j) +h_real[i,j] ==-Q_real[i,j])
            @constraint(opti, second_back_x_ho4(h_imag,dx,i,j) +second_forw_y_ho4(h_imag,dx,i,j) +h_imag[i,j] ==-Q_imag[i,j])
        end
    end      




    # Had to deconstruct the constraints into loops to use the angle finder
    for k in 4:n-3
        # first check angle
        theta = angle_check(k,1,divisions)
        # then add constraint
        #x- real 
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_real,dx,k,1) + sin.(theta)*first_centred_y_ho4(h_real,dx,k,1) == h_imag[k,1]);
        # This  behaves with JuMP
        #x- imag
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_imag,dx,k,1) + sin.(theta)*first_centred_y_ho4(h_imag,dx,k,1) == -h_real[k,1]); 
    end

    for k in 4:n-3
        # first check angle
        theta = angle_check(k,n,divisions)
        # then add constraint
        #x+ real 
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_real,dx,k,n) + sin.(theta)*first_centred_y_ho4(h_real,dx,k,n) == h_imag[k,n]);
        #x+ imag;
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_imag,dx,k,n) + sin.(theta)*first_centred_y_ho4(h_imag,dx,k,n) == -h_real[k,n]); 
    end

    for k in 4:n-3
        # first check angle
        theta = angle_check(n,k,divisions)
        # then add constraint
        #y- real 
        @constraint(opti, cos.(theta)*first_centred_x_ho4(h_real,dx,n,k) + sin.(theta)*first_back_y_ho4(h_real,dx,n,k) == h_imag[n,k]);
        #y- imag
        @constraint(opti, cos.(theta)*first_centred_x_ho4(h_imag,dx,n,k) + sin.(theta)*first_back_y_ho4(h_imag,dx,n,k) == -h_real[n,k]); 
        # swapped sign to see if it fixes
    end

    for k in 4:n-3
        # first check angle
        theta = angle_check(1,k,divisions)
        # then add constraint
        #y+ real 
        @constraint(opti, cos.(theta)*first_centred_x_ho4(h_real,dx,1,k) + sin.(theta)*first_forw_y_ho4(h_real,dx,1,k) == h_imag[1,k]);
        #y+ imag
        @constraint(opti, cos.(theta)*first_centred_x_ho4(h_imag,dx,1,k) + sin.(theta)*first_forw_y_ho4(h_imag,dx,1,k) == -h_real[1,k]);
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
 
    # x+,y-
    for i in 1:3
        theta = angle_check(i,n,divisions)
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_real,dx,i,n) + sin.(theta)*first_forw_y_ho4(h_real,dx,i,n) == h_imag[i,n]);
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_imag,dx,i,n) + sin.(theta)*first_forw_y_ho4(h_imag,dx,i,n) == -h_real[i,n]); 
    end
    for j in n-2:n-1 # be careful to avoid double counting at the corner
        theta = angle_check(1,j,divisions)
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_real,dx,1,j) + sin.(theta)*first_forw_y_ho4(h_real,dx,1,j) == h_imag[1,j]);
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_imag,dx,1,j) + sin.(theta)*first_forw_y_ho4(h_imag,dx,1,j) == -h_real[1,j]); 
    end


    # x+,y+
    for i in n-2:n
        theta = angle_check(i,n,divisions)
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_real,dx,i,n) + sin.(theta)*first_back_y_ho4(h_real,dx,i,n) == h_imag[i,n]);
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_imag,dx,i,n) + sin.(theta)*first_back_y_ho4(h_imag,dx,i,n) == -h_real[i,n]); 
    end
    for i in n-2:n-1
        theta = angle_check(n,i,divisions)
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_real,dx,n,i) + sin.(theta)*first_back_y_ho4(h_real,dx,n,i) == h_imag[n,i]);
        @constraint(opti, cos.(theta)*first_back_x_ho4(h_imag,dx,n,i) + sin.(theta)*first_back_y_ho4(h_imag,dx,n,i) == -h_real[n,i]); 
    end

    # x-,y-
    for i in 1:3
        theta = angle_check(i,1,divisions)
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_real,dx,i,1) + sin.(theta)*first_forw_y_ho4(h_real,dx,i,1) == h_imag[i,1]);
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_imag,dx,i,1) + sin.(theta)*first_forw_y_ho4(h_imag,dx,i,1) == -h_real[i,1]); 
    end
    for i in 2:3
        theta = angle_check(1,i,divisions)
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_real,dx,1,i) + sin.(theta)*first_forw_y_ho4(h_real,dx,1,i) == h_imag[1,i]);
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_imag,dx,1,i) + sin.(theta)*first_forw_y_ho4(h_imag,dx,1,i) == -h_real[1,i]); 
    end

    # x-,y+
    for i in n-2:n-1
        theta = angle_check(i,1,divisions)
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_real,dx,i,1) + sin.(theta)*first_back_y_ho4(h_real,dx,i,1) == h_imag[i,1]);
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_imag,dx,i,1) + sin.(theta)*first_back_y_ho4(h_imag,dx,i,1) == -h_real[i,1]); 
    end
    for i in 1:3
        theta = angle_check(n,i,divisions)
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_real,dx,n,i) + sin.(theta)*first_back_y_ho4(h_real,dx,n,i) == h_imag[n,i]);
        @constraint(opti, cos.(theta)*first_forw_x_ho4(h_imag,dx,n,i) + sin.(theta)*first_back_y_ho4(h_imag,dx,n,i) == -h_real[n,i]); 
    end


    # still has the correct number of constraints

    


    ##################
    ## bounded norm ##
    ##################

    @constraint(opti, dx^2 * sum(sum(Q_real[i,j]^2 
        + Q_imag[i,j]^2 for i in left_index:right_index) for j in left_index:right_index) <= 1);


    # thrust with just the fore-aft
    if solver == false
        # thrust with just the fore-aft
        # @objective(opti, Max, (1/2)*(-(dx/2)*(h_real[left_index,right_index]^2 + h_imag[left_index,right_index]^2 
        #         - h_real[left_index,left_index]^2- h_imag[left_index,left_index]^2 + h_real[right_index,right_index]^2 
        #         + h_imag[right_index,right_index]^2 - h_real[right_index,left_index]^2- h_imag[right_index,left_index]^2)
        #         -dx * sum(h_real[k,right_index]^2 + h_imag[k,right_index]^2 
        #         - h_real[k,left_index]^2 - h_imag[k,left_index]^2 for k in left_index+1:right_index-1)))
        # thrust (that I think is wrong) 
        # @objective(opti, Max, (dx) * sum(sum(Q_real[i,j]*(-h_real[i,j-1] + 
        #             h_real[i,j+1]) + Q_imag[i,j]*(-h_imag[i,j-1] + 
        #             h_imag[i,j+1])  for i in left_index:right_index) 
        #             for j in left_index:right_index));
        @objective(opti, Max,x_thrust_Q(left_index, right_index, left_index, right_index, dx, Q_real, Q_imag,h_real,h_imag))
    else
        @objective(opti, Max, 1);
    end
    optimize!(opti)

    if print_values == true
        ###################################
        ## Calculating various variables ##
        ###################################
        # note these take ages for l_inner == l_outer

        nrm = dx^2 *sum(sum(value.(value.(value.(value.(value.(Q_real)))))[i,j]^2 + value.(Q_imag)[i,j]^2 for i in left_index:right_index) for j in left_index:right_index)
        # pwr = 2*dx^2 * sum(sum(-value.(Q_real)[i,j]*value.(h_imag)[i,j] +
        #                     value.(Q_imag)[i,j]*value.(h_real)[i,j]  for i in left_index:right_index) for j in left_index:right_index)
        # pwr2 = dx * sum((value.(h_imag)[k,left_index])^2 + (value.(h_real)[k,left_index])^2 + 
        #         (value.(h_imag)[k,right_index])^2 + (value.(h_real)[k,right_index])^2 for k in left_index:right_index)
        
        ###############
        ## x related ##
        ###############
        # xforc1 = (dx/4) * sum(sum(value.(Q_real)[i,j]*(-value.(h_real)[i,j-1] + 
        #             value.(h_real)[i,j+1]) + value.(Q_imag)[i,j]*(-value.(h_imag)[i,j-1] + 
        #             value.(h_imag)[i,j+1])  for i in left_index:right_index) for j in left_index:right_index)
        xforc2 = (1/2)*(-(dx/2)*(value.(h_real)[left_index,right_index]^2 + value.(h_imag)[left_index,right_index]^2 
                - value.(h_real)[left_index,left_index]^2- value.(h_imag)[left_index,left_index]^2 + value.(h_real)[right_index,right_index]^2 
                + value.(h_imag)[right_index,right_index]^2 - value.(h_real)[right_index,left_index]^2- value.(h_imag)[right_index,left_index]^2)
                -dx * sum(value.(h_real)[k,right_index]^2 + value.(h_imag)[k,right_index]^2 - value.(h_real)[k,left_index]^2 
                                - value.(h_imag)[k,left_index]^2 for k in left_index+1:right_index-1))
        
        xforc3 = (1/2)*(-(dx/2)*(value.(h_real)[1,right_index]^2 + value.(h_imag)[1,right_index]^2 
                - value.(h_real)[1,left_index]^2- value.(h_imag)[1,left_index]^2 + value.(h_real)[n,right_index]^2 
                + value.(h_imag)[n,right_index]^2 - value.(h_real)[n,left_index]^2- value.(h_imag)[n,left_index]^2)
                -dx * sum(value.(h_real)[k,right_index]^2 + value.(h_imag)[k,right_index]^2 - value.(h_real)[k,left_index]^2 
                                - value.(h_imag)[k,left_index]^2 for k in 2:n-1))

        # Calculating that extra 2d term in the force
        # xforc_extra = 
        xforce_last_term1 = extra_x1(left_index, right_index, left_index, right_index, dx, value.(h_real), value.(h_imag))
        xforce_last_term2 = extra_x2(left_index, right_index, left_index, right_index, dx, value.(h_real), value.(h_imag))
        xforce_Q = x_thrust_Q(left_index, right_index, left_index, right_index, dx, value.(Q_real), value.(Q_imag),value.(h_real),value.(h_imag))


        ###############
        ## y related ##
        ###############
        # yforce_Q = y_thrust_Q(right_index, left_index, left_index, right_index, dx, value.(Q_real), value.(Q_imag),value.(h_real),value.(h_imag))
        # yforc_foreaft = y_thrust_fore_aft(right_index, left_index, left_index, right_index, dx, value.(h_real), value.(h_imag))
        # yforce_last_term1 = extra_y1(right_index, left_index, left_index, right_index, dx, value.(h_real), value.(h_imag))
        # yforce_last_term2 = extra_y2(right_index, left_index, left_index, right_index, dx, value.(h_real), value.(h_imag))

        
        print(" ","\n left = ", left, 
        "\n left_index = ", left_index,
        "\n right = ", right,
        "\n right index = ", right_index,"\n",
        "nrm check = ", nrm,"\n")


        print("\n xforce (Q) = ", xforce_Q, 
        "\n xforce (fore-aft) = ", xforc2, 
        "\n xforce (fore-aft whole y) = ", xforc3,
        "\n xforce (other term 1) = ", xforce_last_term1,
        "\n xforce (other term 2) = ", xforce_last_term2,
        "\n adding xforces = ", xforce_Q + xforce_last_term1+xforce_last_term2#matches eq 47 in code_explnation2
        )# ,"\n yforce (Q) = ", yforce_Q,
        # "\n yforce (fore-aft term) = ", yforc_foreaft/2,
        # "\n yforce (other term 1) = ", yforce_last_term1,
        # "\n yforce (other term 2) = ", yforce_last_term2)
    end

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

    return q_real, q_imag, H_real, H_imag, dx, left_index, right_index
end
