# modify sup_ellipse


using JuMP
import Ipopt
using Plots
include(relpath("functions/integral_functions.jl","shape_in_wdp"))
include(relpath("functions/finite_functions.jl","shape_in_wdp"))

# Functions for the matrices after modification to be for direct input
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

# Crank-nicolson functions - not used yet
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


# now build the code 
# dimensions:
# the outer region is 15pi/2 as per usual
# inner region 3pi/2 will be hard coded for starters 
# len = 12
function sup_param(body_len_x, body_len_y, v, divs=100, mode=1, extend_domain=false, sil=false,nrm=1,grad_bound=false,solver=false)
    # given the body length, we will calculate the sliced domain so that the 
    # wave reaches the wall at least 2 times the body length behind the body
    ymax = round(2*body_len_x/sqrt(v^2-1) + body_len_y/2,digits=2)
    ymin = round(-2*body_len_x/sqrt(v^2-1) - body_len_y/2,digits=2) 

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
    opti = Model(Ipopt.Optimizer)
    if sil==true
        set_silent(opti)
    end
    # define the variables 
    @variable(opti, h_real[1:n,1:end_t], start=0.1);
    @variable(opti, h_imag[1:n,1:end_t], start=0.1);
    @variable(opti, u_real[1:n,1:end_t], start=0.1);
    @variable(opti, u_imag[1:n,1:end_t], start=0.1);
    if solver == false
        @variable(opti, Q_real[1:n,1:end_t], start=0.1);
        @variable(opti, Q_imag[1:n,1:end_t], start=0.1);
        @variable(opti, A_real,start=0.003);
        @variable(opti, A_imag,start=0.008);
        @variable(opti, B_real,start=0.08258);
        @variable(opti, B_imag,start=0.09822);
        @variable(opti, C_real,start=0.00706);
        @variable(opti, C_imag,start=0.01248);
        @variable(opti, D_real,start=0.08208);
        @variable(opti, D_imag,start=0.09806);
        if mode ==1
            @variable(opti, -1<=k_num<=1,start=0.88204);
        elseif mode == 3
            @variable(opti, 1<=k_num,start=sqrt(v/(v^2-1)));
        elseif mode == 4
            @variable(opti, 1<=k_num,start=1);
        end
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
                    Q_real[yi,ti] = 0.5*cos.(t);
                    Q_imag[yi,ti] = -0.5*sin.(t); 
        # x=-t does not flip sign because it is already baked into t = rhs_t 
                else
                    if mode ==1
                    @constraint(opti, Q_real[yi,ti] == 
                A_real*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            -A_imag*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +B_real*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            -B_imag*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            +C_real*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            -C_imag*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +D_real*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            -D_imag*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y));
                @constraint(opti, Q_imag[yi,ti] == 
                A_imag*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +A_real*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +B_imag*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            +B_real*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            +C_imag*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +C_real*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t +sqrt(1-k_num^2)*y)
            +D_imag*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y)
            +D_real*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t -sqrt(1-k_num^2)*y));
            #         #2a
                    elseif mode == 3
                    @constraint(opti, Q_real[yi,ti] == 
                A_real*cos.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            -A_imag*sin.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +B_real*cos.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            -B_imag*sin.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +C_real*cos.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            -C_imag*sin.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +D_real*cos.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            -D_imag*sin.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1)))));
                @constraint(opti, Q_imag[yi,ti] == 
                A_imag*cos.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +A_real*sin.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +B_imag*cos.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +B_real*sin.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y -t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +C_imag*cos.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +C_real*sin.((v/(1-v^2))*t)*exp(+sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +D_imag*cos.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1))))
            +D_real*sin.((v/(1-v^2))*t)*exp(-sqrt(k_num^2-1)*y +t*sqrt(-(v^2/(1-v^2)^2)+((k_num^2)/(v^2-1)))));
            #         #2b
                    elseif mode == 4
                    @constraint(opti, Q_real[yi,ti] == 
                A_real*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            -A_imag*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +B_real*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            -B_imag*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            +C_real*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            -C_imag*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +D_real*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            -D_imag*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y));
                @constraint(opti, Q_imag[yi,ti] == 
                A_imag*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +A_real*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +B_imag*cos.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            +B_real*sin.((v/(1-v^2) - sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            +C_imag*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +C_real*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(+sqrt(k_num^2-1)*y)
            +D_imag*cos.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y)
            +D_real*sin.((v/(1-v^2) + sqrt((v^2/(1-v^2)^2)+((k_num^2)/(1-v^2))))*t)*exp(-sqrt(k_num^2-1)*y));
                    end
                end
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

    # need bounded norm constraint as well as the objective 
    # fix integrals first. 
    # bounded norm constraint
    if solver == false
        @constraint(opti, dt*dy * sum(sum(Q_real[i,j]^2 + Q_imag[i,j]^2 
                    for i in lhs_index-1:rhs_index+1) for j in 1:body_t) <= nrm)#*pi*(body_len_x/2)*(body_len_y/2));
        # @constraint(opti, sup_q_pwr_sign_flip(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag,v)<=1);

        @objective(opti, Max, -sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag));
        # also tried start index of 2 and it still blows up.
        # constrain thrust to be less than power
        # @constraint(opti, sup_q_pwr_sign_flip(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag,v) >= sup_x_thrust_Q(lhs_index, rhs_index, 1, body_t, dt, dy,Q_real, Q_imag,h_real,h_imag));
    else
        @objective(opti, Max, 1);
    end


    optimize!(opti)

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

        values_to_output = [value.(A_real),value.(A_imag),value.(B_real),value.(B_imag),value.(C_real),value.(C_imag),
                        value.(D_real),value.(D_imag),value.(k_num)]


        # return different matrices
        return qr_new, qi_new, hr_new, hi_new, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t,pwr,values_to_output
    end
    values_to_output = [value.(A_real),value.(A_imag),value.(B_real),value.(B_imag),value.(C_real),value.(C_imag),
                        value.(D_real),value.(D_imag),value.(k_num)]


    return qr, qi, hr, hi, tvals, yvals, thrust, dt, dy,lhs_index,rhs_index,body_t,pwr,values_to_output, is_optimal, is_global
end
