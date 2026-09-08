# working supercritical optimisation 
# This Crank-Nicolson scheme does a coordinate transformation x = -t to march 
# forward in time rather than backwards in x
using JuMP
import Ipopt
using Plots
include(relpath("functions/integral_functions.jl","shape-in-wdp"))
include(relpath("functions/finite_functions.jl","shape-in-wdp"))


"""
    Defines the 2d regularised Dirac delta distribution 
    for prescribed source.
"""
function delta_func(a,y,x)
    return (1/(2*pi*(a^2)))*exp(-(x^2+y^2)/(2*a^2))
end
# functions to fill in the crank-nicolson lhs and rhs as equality constraints
function one(hr, ur, yi, ti, dt)
    return hr[yi,ti] - dt*ur[yi,ti]/2
end 

function two(hi, ui, yi, ti, dt)
    return hi[yi,ti] - dt*ui[yi,ti]/2
end

function three(hr, ur, ui, yi, ti, dy, dt, v)
    return (1-v^2)*ur[yi,ti]/dt - v*ui[yi,ti] + hr[yi,ti]/2 + (hr[yi+1,ti] - 2*hr[yi,ti] + hr[yi-1,ti])/(2*dy^2)
end

function four(hi, ur, ui, yi, ti, dy, dt, v)
    return (1-v^2)*ui[yi,ti]/dt + v*ur[yi,ti] + hi[yi,ti]/2 + (hi[yi+1,ti] - 2*hi[yi,ti] + hi[yi-1,ti])/(2*dy^2)
end

# Also will fill in the right hand side with functions
function rhs_one(yi, ti, hr, ur, dt)
    return hr[yi,ti-1] +ur[yi,ti-1]*dt/2 
end

function rhs_two(yi, ti, hi, ui, dt)
    return hi[yi,ti-1] +ui[yi,ti-1]*dt/2 
end

function rhs_three(yi, ti, hr, ui, ur, Qr, dt, dy, v)
    res = (1-v^2)/dt*ur[yi,ti-1]
    res += v*ui[yi,ti-1]
    res += -hr[yi,ti-1]/2
    res += -Qr[yi,ti]/2
    res += -Qr[yi,ti-1]/2
    res += -hr[yi-1,ti-1]/(2*dy^2)
    res += hr[yi,ti-1]/(dy^2)
    res += -hr[yi+1,ti-1]/(2*dy^2)
    return res
end

function rhs_four(yi, ti, hi, ui, ur, Qi, dt, dy, v)
    res = (1-v^2)*ui[yi,ti-1]/dt
    res += -v*ur[yi,ti-1]
    res += -hi[yi,ti-1]/2
    res += -Qi[yi,ti]/2
    res += -Qi[yi,ti-1]/2
    res += -hi[yi-1,ti-1]/(2*dy^2)
    res += hi[yi,ti-1]/(dy^2)
    res += -hi[yi+1,ti-1]/(2*dy^2)
    return res
end

# Crank-Nicolson functions used to extend the domain

function cn_one(A, j, dt, n)
    A[j,j] = 1
    A[j,j+2n] = -dt/2
end 

function cn_two(A, j, dt, n)
    A[j+n,j+n] = 1
    A[j+n,j+3n] = -dt/2
end

function cn_three(A, j, dt, dy, v, n)
    A[j+2n,j+2n] = (1-v^2)/dt
    A[j+2n,j+3n] = -v
    A[j+2n,j] = 1/2 - 1/dy^2
    A[j+2n,j-1] = 1/(2*dy^2)
    A[j+2n,j+1] = 1/(2*dy^2)
end

function cn_four(A, j, dt, dy, v, n)
    A[j+3n,j+3n] = (1-v^2)/dt
    A[j+3n,j+2n] = v
    A[j+3n,j+n] = 1/2 - 1/dy^2
    A[j+3n,j-1+n] = 1/(2*dy^2)
    A[j+3n,j+1+n] = 1/(2*dy^2)
end

# Also will fill in the right hand side with functions
function cn_rhs_one(yi, ti, rhs, hr, ur, dt)
    rhs[yi] = hr[yi,ti-1] +ur[yi,ti-1]*dt/2 
end

function cn_rhs_two(yi, ti, rhs, hi, ui, dt, n)
    rhs[yi+n] = hi[yi,ti-1] +ui[yi,ti-1]*dt/2 
end

function cn_rhs_three(yi, ti, rhs, hr, ui, ur, Qr, dt, dy, v, n)
    res = (1-v^2)/dt*ur[yi,ti-1]
    res += v*ui[yi,ti-1]
    res += -hr[yi,ti-1]/2
    res += -Qr[yi,ti]/2
    res += -Qr[yi,ti-1]/2
    res += -hr[yi-1,ti-1]/(2*dy^2)
    res += hr[yi,ti-1]/(dy^2)
    res += -hr[yi+1,ti-1]/(2*dy^2)
    rhs[yi+2n] = res
end

function cn_rhs_four(yi, ti, rhs, hi, ui, ur, Qi, dt, dy, v, n)
    res = (1-v^2)*ui[yi,ti-1]/dt
    res += -v*ur[yi,ti-1]
    res += -hi[yi,ti-1]/2
    res += -Qi[yi,ti]/2
    res += -Qi[yi,ti-1]/2
    res += -hi[yi-1,ti-1]/(2*dy^2)
    res += hi[yi,ti-1]/(dy^2)
    res += -hi[yi+1,ti-1]/(2*dy^2)
    rhs[yi+3n] = res
end

# optimisation function
function sup_opti(body_len_x, body_len_y, v, divs=100, extend_domain=false, sil=false,solver=false)
    # given the body length, we will calculate the sliced domain so that the 
    # wave reaches the wall at least 4 times the body length behind the body
    ymax = round(2*body_len_x/sqrt(v^2-1) + body_len_y/2,digits=2)
    ymin = round(-2*body_len_x/sqrt(v^2-1) - body_len_y/2,digits=2) 
    #### This 2 times the body length can be tuned if you are extending the 
    #### domain to 4 times which makes the body narrower in points 
    #### for the same number of divisions.
    n = divs + 1
    dy = (ymax-ymin)/divs
    dt =  0.1*dy # make it smaller

    # finding the body region 
    middle = Int64(divs/2 +1)
    body_y = Int64(floor(body_len_y/(2*dy)))
    body_t = Int64(floor(body_len_x/(dt))) # body length is equal in both directions 
    lhs_index = middle-body_y
    rhs_index = middle+body_y
    # time length

    # if extend_domain == false
    end_t = body_t + 5
    # else
    #     # end_t = 4*body_t - 5 
    # end
    # only need to be outside of the body, no worries about visualisation 
    tvals = [];
    yvals = [];
    for i in 1:n
        append!(yvals, ymin+(i-1)*dy)
    end
    for i in 1:end_t
        append!(tvals, (i-1)*dt)
    end


    # We will begin with a square region of length 3pi/2 
    # solver = true # engage the solver rather than optimisation 
    opti = Model(Ipopt.Optimizer)
    if sil==true
        set_silent(opti)
    end
    # set_attribute(opti, "max_cpu_time", 30.0)
    # define the variables 
    @variable(opti, h_real[1:n,1:end_t], start=0.1);
    @variable(opti, h_imag[1:n,1:end_t], start=0.1);
    @variable(opti, u_real[1:n,1:end_t], start=0.1);
    @variable(opti, u_imag[1:n,1:end_t], start=0.1);
    if solver == false
        @variable(opti, Q_real[1:n,1:end_t], start=0.1);
        @variable(opti, Q_imag[1:n,1:end_t], start=0.1);
    else
        Q_real = zeros(n,end_t);
        Q_imag = zeros(n,end_t);
    end
    # impose the source only being defined in a certain region
    ###########################################
    ## Forcing it to be a elliptical shape ####
    # ======================================= #
    
    a = body_len_x/2
    b = body_len_y/2
    rhs_t = a + dt 
    for ti in 1:end_t
        for yi in 1:n
            t = rhs_t - (ti-1)*dt
            y = ymax - (yi-1)*dy
            if (t/a)^2 + (y/b)^2 > 1
                if solver == false
                    @constraint(opti, Q_real[yi,ti] == 0);
                    @constraint(opti, Q_imag[yi,ti] == 0);
                end
            else
                if solver  == true
                    Q_real[yi,ti] = 0.5*cos.(t/(1+v));
                    Q_imag[yi,ti] = -0.5*sin.(t/(1+v)); 
        # x=-t does not flip sign because it is already baked into t = rhs_t ...
                end
            end
        end
    end

    # Crank-Nicolson in the bulk
    for ti in 2:end_t
        for yi in 1:n
            if yi == 1 || yi == n
                @constraint(opti, h_real[yi,ti] == 0);
                @constraint(opti, h_imag[yi,ti] == 0);
                @constraint(opti, u_real[yi,ti] == 0);
                @constraint(opti, u_imag[yi,ti] == 0);
            else
                @constraint(opti, one(h_real, u_real, yi, ti, dt) == rhs_one(yi, ti, h_real, u_real, dt))
                @constraint(opti, two(h_imag, u_imag, yi, ti, dt) == rhs_two(yi, ti, h_imag, u_imag, dt))
                @constraint(opti, three(h_real, u_real, u_imag, yi, ti, dy, dt, v) == rhs_three(yi, ti, h_real, u_imag, u_real, Q_real, dt, dy, v))
                @constraint(opti, four(h_imag, u_real, u_imag, yi, ti, dy, dt, v) == rhs_four(yi, ti, h_imag, u_imag, u_real, Q_imag, dt, dy, v))
            end
        end
    end
    # need to constrain the first index to be zero (initial condition)
    for i in 1:n
        @constraint(opti, h_real[i,1] == 0);
        @constraint(opti, h_imag[i,1] == 0);
        @constraint(opti, u_real[i,1] == 0);
        @constraint(opti, u_imag[i,1] == 0);
    end
    # need bounded norm constraint as well as the objective 
    if solver == false
        @constraint(opti, dt*dy * sum(sum(Q_real[i,j]^2 + Q_imag[i,j]^2 for i in lhs_index-1:rhs_index+1) for j in 1:body_t) <= 1)#*pi*(body_len_x/2)*(body_len_y/2));
        @objective(opti, Max, -sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag)); 
        # sign change in thrust due to the coordinate transformation x = -t
        # constrain thrust to be less than power in case things break
        @constraint(opti, sup_q_pwr_sign_flip(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag,v) >= -sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag));
    else
        @objective(opti, Max, 1);
    end

    optimize!(opti)

    # exporting values
    hr = value.(h_real)
    hi = value.(h_imag)
    qr = value.(Q_real)
    qi = value.(Q_imag)
    ur = value.(u_real)
    ui = value.(u_imag)

    thrust = -sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy, qr, qi,hr,hi)
    pwr = sup_q_pwr_sign_flip(lhs_index, rhs_index, 1, body_t, dt, dy,qr, qi,hr,hi,v)
    is_optimal = is_solved_and_feasible(opti)
    is_global =  is_solved_and_feasible(opti; allow_local=false)
    # put in an if statement for the extended domain with C-N
    if extend_domain == true
        end_t = 4*body_t 
        # extend the matrices by concatenation 
        hr_new = hcat(hr, zeros(Float64, size(hr, 1), 3*body_t-5))
        hi_new = hcat(hi, zeros(Float64, size(hi, 1), 3*body_t-5))
        qr_new = hcat(qr, zeros(Float64, size(qr, 1), 3*body_t-5))
        qi_new = hcat(qi, zeros(Float64, size(qi, 1), 3*body_t-5))
        ur_new = hcat(ur, zeros(Float64, size(ur, 1), 3*body_t-5))
        ui_new = hcat(ui, zeros(Float64, size(ui, 1), 3*body_t-5))
        # set up the Crank-Nicolson solver
        A = zeros(4n,4n);

        for yi in 2:n-1
            cn_one(A, yi, dt, n)
            cn_two(A, yi, dt, n)
            cn_three(A, yi, dt, dy, v, n)
            cn_four(A, yi, dt, dy, v, n)
        end
        # boundary conditions
        A[1,1] = 1;
        A[n,n] = 1;
        A[n+1,n+1] = 1;
        A[2n,2n] = 1;
        A[2n+1,2n+1] = 1;
        A[3n,3n] = 1;
        A[3n+1,3n+1] = 1;
        A[4n,4n] = 1;

        invA = inv(A);

        for ti in body_t+6:end_t
            rhs = zeros(4n,1)
            for yi in 2:n-1
                cn_rhs_one(yi, ti, rhs, hr_new, ur_new, dt)
                cn_rhs_two(yi, ti, rhs, hi_new, ui_new, dt, n)
                cn_rhs_three(yi,ti,rhs,hr_new,ui_new,ur_new,qr_new,dt,dy,v,n)
                cn_rhs_four(yi,ti,rhs,hi_new,ui_new,ur_new,qi_new,dt,dy,v,n)
            end
            x = invA*rhs
            # fill in the matrices
            hr_new[:,ti] = x[1:n]
            hi_new[:,ti] = x[n+1:2n]
            ur_new[:,ti] = x[2n+1:3n]
            ui_new[:,ti] = x[3n+1:4n]
        end
        # redo the t value array 
        tvals = []
        for i in 1:end_t
            append!(tvals, (i-1)*dt)
        end
        # return different matrices
        return qr_new, qi_new, hr_new, hi_new, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t,pwr,ur_new,ui_new
    end
    return qr, qi, hr, hi, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t,pwr,is_optimal,is_global, ur, ui
end


function prescribed_source(body_len_x, body_len_y, v, divs=100, extend_domain=false, sil=false)
    # body_len = 4
    # v=2
    # given the body length, we will calculate the sliced domain so that the 
    # wave reaches the wall at least 2 times the body length behind the body
    ymax = round(4*body_len_x/sqrt(v^2-1) + body_len_y/2,digits=2)
    ymin = round(-4*body_len_x/sqrt(v^2-1) - body_len_y/2,digits=2)
    
 
    n = divs + 1
    dy = (ymax-ymin)/divs
    dt =  0.1*dy # make it smaller

    # finding the body region 
    middle = Int64(divs/2 +1)
    body_y = Int64(floor(body_len_y/(2*dy)))
    body_t = Int64(floor(body_len_x/(dt))) # body length is equal in both directions 
    lhs_index = middle-body_y
    rhs_index = middle+body_y
    # time length
    end_t = body_t + 5
    tvals = [];
    yvals = [];
    for i in 1:n
        append!(yvals, ymin+(i-1)*dy)
    end
    for i in 1:end_t
        append!(tvals, (i-1)*dt)
    end


    # We will begin with a square region of length 3pi/2 
    solver = false # engage the solver rather than optimisation 
    opti = Model(Ipopt.Optimizer)
    if sil==true
        set_silent(opti)
    end
    # set_attribute(opti, "max_cpu_time", 30.0)
    # define the variables 
    @variable(opti, h_real[1:n,1:end_t], start=0.1);
    @variable(opti, h_imag[1:n,1:end_t], start=0.1);
    @variable(opti, u_real[1:n,1:end_t], start=0.1);
    @variable(opti, u_imag[1:n,1:end_t], start=0.1);
    Q_real = zeros(n,end_t);
    Q_imag = zeros(n,end_t);

    # impose the source only being defined in a certain region
    ###########################################
    ## Forcing it to be a elliptical shape ####
    # ======================================= #
    
    a = body_len_x/2
    b = body_len_y/2
    rhs_t = a + dt 
    for ti in 1:end_t
        for yi in 1:n
            t = rhs_t - (ti-1)*dt
            y = ymax - (yi-1)*dy
            if (t/a)^2 + (y/b)^2 <= 1
                Q_real[yi,ti] = delta_func(dy,y,t)
            end
        end
    end
    

    # build a matrix for the right hand side, then do some multiplications 
    # If this fails, I will apply them manually 
    for ti in 2:end_t
        for yi in 1:n
            if yi == 1 || yi == n
                @constraint(opti, h_real[yi,ti] == 0);
                @constraint(opti, h_imag[yi,ti] == 0);
                @constraint(opti, u_real[yi,ti] == 0);
                @constraint(opti, u_imag[yi,ti] == 0);
            else
                @constraint(opti, one(h_real, u_real, yi, ti, dt) == rhs_one(yi, ti, h_real, u_real, dt))
                @constraint(opti, two(h_imag, u_imag, yi, ti, dt) == rhs_two(yi, ti, h_imag, u_imag, dt))
                @constraint(opti, three(h_real, u_real, u_imag, yi, ti, dy, dt, v) == rhs_three(yi, ti, h_real, u_imag, u_real, Q_real, dt, dy, v))
                @constraint(opti, four(h_imag, u_real, u_imag, yi, ti, dy, dt, v) == rhs_four(yi, ti, h_imag, u_imag, u_real, Q_imag, dt, dy, v))
            end
        end
        # print(ti,"\t")
    end
    # need to constrain the first index to be zero
    for i in 1:n
        @constraint(opti, h_real[i,1] == 0);
        @constraint(opti, h_imag[i,1] == 0);
        @constraint(opti, u_real[i,1] == 0);
        @constraint(opti, u_imag[i,1] == 0);
    end

    # this should be all the constraints then 

    # trivial optimisation
    @objective(opti, Max, 1);



    optimize!(opti)

    hr = value.(h_real)
    hi = value.(h_imag)
    qr = value.(Q_real)
    qi = value.(Q_imag)
    ur = value.(u_real)
    ui = value.(u_imag)

    thrust = sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy, qr, qi,hr,hi)

    # put in an if statement for the extended domain with C-N
    if extend_domain == true
        end_t = 4*body_t 
        # extend the matrices by concatenation 
        hr_new = hcat(hr, zeros(Float64, size(hr, 1), 3*body_t-5))
        hi_new = hcat(hi, zeros(Float64, size(hi, 1), 3*body_t-5))
        qr_new = hcat(qr, zeros(Float64, size(qr, 1), 3*body_t-5))
        qi_new = hcat(qi, zeros(Float64, size(qi, 1), 3*body_t-5))
        ur_new = hcat(ur, zeros(Float64, size(ur, 1), 3*body_t-5))
        ui_new = hcat(ui, zeros(Float64, size(ui, 1), 3*body_t-5))
        # set up the Crank-Nicolson solver
        A = zeros(4n,4n);

        for yi in 2:n-1
            cn_one(A, yi, dt, n)
            cn_two(A, yi, dt, n)
            cn_three(A, yi, dt, dy, v, n)
            cn_four(A, yi, dt, dy, v, n)
        end
        # boundary conditions
        A[1,1] = 1;
        A[n,n] = 1;
        A[n+1,n+1] = 1;
        A[2n,2n] = 1;
        A[2n+1,2n+1] = 1;
        A[3n,3n] = 1;
        A[3n+1,3n+1] = 1;
        A[4n,4n] = 1;

        invA = inv(A);

        for ti in body_t+6:end_t
            rhs = zeros(4n,1)
            for yi in 2:n-1
                cn_rhs_one(yi, ti, rhs, hr_new, ur_new, dt)
                cn_rhs_two(yi, ti, rhs, hi_new, ui_new, dt, n)
                cn_rhs_three(yi,ti,rhs,hr_new,ui_new,ur_new,qr_new,dt,dy,v,n)
                cn_rhs_four(yi,ti,rhs,hi_new,ui_new,ur_new,qi_new,dt,dy,v,n)
            end
            x = invA*rhs
            # fill in the matrices
            hr_new[:,ti] = x[1:n]
            hi_new[:,ti] = x[n+1:2n]
            ur_new[:,ti] = x[2n+1:3n]
            ui_new[:,ti] = x[3n+1:4n]
        end
        # redo the t value array 
        tvals = []
        for i in 1:end_t
            append!(tvals, (i-1)*dt)
        end
        # return different matrices
        return qr_new, qi_new, hr_new, hi_new, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t
    end
    return qr, qi, hr, hi, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t
end
