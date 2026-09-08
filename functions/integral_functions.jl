# storage for all the integral functions
include("finite_functions.jl")
# x integrals


function extra_x1(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   #=
       Takes a matrix A in its real part and imaginary part and does the 
       integral I^x_1 in the notes in the "2D Tasks/Studying force terms" page
       in the square within the bounds of yminus, yplus, xminus, xplus inputs
   =#
   ans = - (dx/4) * (Areal[yminus,xplus]*first_back_x(Aimag, dx, yminus,xplus)
               - Aimag[yminus,xplus]*first_back_x(Areal, dx, yminus,xplus)
               + Areal[yplus,xplus]*first_back_x(Aimag, dx, yplus,xplus)
               - Aimag[yplus,xplus]*first_back_x(Areal, dx, yplus,xplus))
   ans += - (dx/4) * (Areal[yminus,xminus]*first_forw_x(Aimag, dx, yminus,xminus)
                   - Aimag[yminus,xminus]*first_forw_x(Areal, dx, yminus,xminus)
                   + Areal[yplus,xminus]*first_forw_x(Aimag, dx, yplus,xminus)
                   - Aimag[yplus,xminus]*first_forw_x(Areal, dx, yplus,xminus))
   ans += - (dx/2) * sum(Areal[yminus,k]*first_centred_x(Aimag, dx, yminus,k)
                   - Aimag[yminus,k]*first_centred_x(Areal, dx, yminus,k)
                   + Areal[yplus,k]*first_centred_x(Aimag, dx, yplus,k)
                   - Aimag[yplus,k]*first_centred_x(Areal, dx, yplus,k) 
                   for k in xminus+1:xplus-1)
   return ans
end


"""
   Made a modified version to work with when we are adjusting y
"""
function extra_x1_mod(yminus, yplus, xminus, xplus, dx, Areal, Aimag,coef=1)
   #=
       Takes a matrix A in its real part and imaginary part and does the 
       integral I^x_1 in the notes in the "2D Tasks/Studying force terms" page
       in the square within the bounds of yminus, yplus, xminus, xplus inputs
   =#
   ans = - coef*(dx/4) * (Areal[yminus,xplus]*first_back_x(Aimag, dx, yminus,xplus)
               - Aimag[yminus,xplus]*first_back_x(Areal, dx, yminus,xplus)
               + Areal[yplus,xplus]*first_back_x(Aimag, dx, yplus,xplus)
               - Aimag[yplus,xplus]*first_back_x(Areal, dx, yplus,xplus))
   ans += - coef*(dx/4) * (Areal[yminus,xminus]*first_forw_x(Aimag, dx, yminus,xminus)
                   - Aimag[yminus,xminus]*first_forw_x(Areal, dx, yminus,xminus)
                   + Areal[yplus,xminus]*first_forw_x(Aimag, dx, yplus,xminus)
                   - Aimag[yplus,xminus]*first_forw_x(Areal, dx, yplus,xminus))
   ans += - coef*(dx/2) * sum(Areal[yminus,k]*first_centred_x(Aimag, dx, yminus,k)
                   - Aimag[yminus,k]*first_centred_x(Areal, dx, yminus,k)
                   + Areal[yplus,k]*first_centred_x(Aimag, dx, yplus,k)
                   - Aimag[yplus,k]*first_centred_x(Areal, dx, yplus,k) 
                   for k in xminus+1:xplus-1)
   return ans
end

function extra_x2(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   #=
       Takes a matrix A in its real part and imaginary part and does the 
       integral I^y_2 in the notes in the square within the bounds of 
       yminus, yplus, xminus, xplus inputs

   =#
   ans2 = - (dx/8) * (first_back_y(Areal, dx, yplus, xplus)^2 
                   + first_back_y(Aimag, dx, yplus, xplus)^2
                   - first_back_y(Areal, dx, yplus, xminus)^2 
                   - first_back_y(Aimag, dx, yplus, xminus)^2
                   + first_forw_y(Areal, dx, yminus, xplus)^2 
                   + first_forw_y(Aimag, dx, yminus, xplus)^2
                   - first_forw_y(Areal, dx, yminus, xminus)^2 
                   - first_forw_y(Aimag, dx, yminus, xminus)^2 )
   ans2 += - (dx/4) * sum(first_centred_y(Areal, dx, k, xplus)^2 
                   + first_centred_y(Aimag, dx, k, xplus)^2
                   - first_centred_y(Areal, dx, k, xminus)^2 
                   - first_centred_y(Aimag, dx, k, xminus)^2
                   for k in yminus+1:yplus-1)
   return ans2
end


"""
Calculates the thrust term with the Q and h_x
"""
function x_thrust_Q(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag)
   # corners
   ans = (1/8)*(dx^2)*(qreal[yplus,xminus]*first_forw_x(hreal,dx,yplus,xminus) +
                        qimag[yplus,xminus]*first_forw_x(himag,dx,yplus,xminus) +
                        qreal[yplus,xplus]*first_back_x(hreal,dx,yplus,xplus) +
                        qimag[yplus,xplus]*first_back_x(himag,dx,yplus,xplus) +
                        qreal[yminus,xminus]*first_forw_x(hreal,dx,yminus,xminus)+
                        qimag[yminus,xminus]*first_forw_x(himag,dx,yminus,xminus)+
                        qreal[yminus,xplus]*first_back_x(hreal,dx,yminus,xplus) +
                        qimag[yminus,xplus]*first_back_x(himag,dx,yminus,xplus) )
   # edges
   ans += (1/4)*(dx^2)*(sum( qreal[k,xminus]*first_forw_x(hreal,dx,k,xminus)
                           + qimag[k,xminus]*first_forw_x(himag,dx,k,xminus)
                           + qreal[k,xplus]*first_back_x(hreal,dx,k,xplus)
                           + qimag[k,xplus]*first_back_x(himag,dx,k,xplus)
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*first_centred_x(hreal,dx,yplus,k)
                           + qimag[yplus,k]*first_centred_x(himag,dx,yplus,k)
                           + qreal[yminus,k]*first_centred_x(hreal,dx,yminus,k)
                           + qimag[yminus,k]*first_centred_x(himag,dx,yminus,k)
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx^2)*(sum( sum( qreal[k,m]*first_centred_x(hreal,dx,k,m)
                                 + qimag[k,m]*first_centred_x(himag,dx,k,m)
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )
   return ans
end 

# y

# function y_thrust_Q(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag)
#     #= 
#         Explanation
#     =#
#     return ((dx) * sum(sum(qreal[i,j]*(-hreal[i-1,j] + hreal[i+1,j]) + qimag[i,j]*(-himag[i-1,j] + himag[i+1,j])  for i in yminus:yplus) for j in xminus:xplus) )
# end
# derivatives might be backwards because the indices are inverted in matrices
# This might be the wrong approach

"""
Calculates the thrust term in y with the Q and h_y
"""
function y_thrust_Q(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag)
   # corners
   ans = (1/8)*(dx^2)*(qreal[yplus,xminus]*first_back_y(hreal,dx,yplus,xminus) +
                        qimag[yplus,xminus]*first_back_y(himag,dx,yplus,xminus) +
                        qreal[yplus,xplus]*first_back_y(hreal,dx,yplus,xplus) +
                        qimag[yplus,xplus]*first_back_y(himag,dx,yplus,xplus) +
                        qreal[yminus,xminus]*first_forw_y(hreal,dx,yminus,xminus)+
                        qimag[yminus,xminus]*first_forw_y(himag,dx,yminus,xminus)+
                        qreal[yminus,xplus]*first_forw_y(hreal,dx,yminus,xplus) +
                        qimag[yminus,xplus]*first_forw_y(himag,dx,yminus,xplus) )
   # edges
   ans += (1/4)*(dx^2)*(sum( qreal[k,xminus]*first_centred_y(hreal,dx,k,xminus)
                           + qimag[k,xminus]*first_centred_y(himag,dx,k,xminus)
                           + qreal[k,xplus]*first_centred_y(hreal,dx,k,xplus)
                           + qimag[k,xplus]*first_centred_y(himag,dx,k,xplus)
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*first_back_y(hreal,dx,yplus,k)
                           + qimag[yplus,k]*first_back_y(himag,dx,yplus,k)
                           + qreal[yminus,k]*first_forw_y(hreal,dx,yminus,k)
                           + qimag[yminus,k]*first_forw_y(himag,dx,yminus,k)
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx^2)*(sum( sum( qreal[k,m]*first_centred_y(hreal,dx,k,m)
                                 + qimag[k,m]*first_centred_y(himag,dx,k,m)
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )
   return ans
end 


function y_thrust_fore_aft(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
    #= 
        Explanation
    =#
    ###############################################
    ## I think there is a factor of 1/2 missing ###
    ###############################################
    return (-(dx/2)*(Areal[yplus,xplus]^2 + Aimag[yplus,xplus]^2 
            - Areal[yminus,xplus]^2- Aimag[yminus,xplus]^2 
            + Areal[yplus,xminus]^2 + Aimag[yplus,xminus]^2 
            - Areal[yminus,xminus]^2- Aimag[yminus,xminus]^2) -
            dx * sum(Areal[yplus,k]^2 + Aimag[yplus,k]^2 
            - Areal[yminus,k]^2 - Aimag[yminus,k]^2 for k in yplus+1:yminus-1)) 
end

function extra_y1(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
    #=
        Takes a matrix A in its real part and imaginary part and does the 
        integral I^y_1 in the notes in the "2D Tasks/Studying force terms" page
        in the square within the bounds of yminus, yplus, xminus, xplus inputs
    =#
    ans = - (dx/4) * (Areal[yplus,xplus]*first_forw_y(Aimag, dx, yplus,xplus)
                - Aimag[yplus,xplus]*first_forw_y(Areal, dx, yplus,xplus)
                + Areal[yplus,xminus]*first_forw_y(Aimag, dx, yplus,xminus)
                - Aimag[yplus,xminus]*first_forw_y(Areal, dx, yplus,xminus))
    ans += - (dx/4) * (Areal[yminus,xplus]*first_back_y(Aimag, dx, yminus,xplus)
                    - Aimag[yminus,xplus]*first_back_y(Areal, dx, yminus,xplus)
                    + Areal[yminus,xminus]*first_back_y(Aimag, dx, yminus,xminus)
                    - Aimag[yminus,xminus]*first_back_y(Areal, dx, yminus,xminus))
    ans += - (dx/2) * sum(Areal[k,xplus]*first_centred_y(Aimag, dx, k,xplus)
                    - Aimag[k,xplus]*first_centred_y(Areal, dx, k,xplus)
                    + Areal[k,xminus]*first_centred_y(Aimag, dx, k,xminus)
                    - Aimag[k,xminus]*first_centred_y(Areal, dx, k,xminus) 
                    for k in yplus+1:yminus-1)
    return ans
end

function extra_y2(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
    #=
        Takes a matrix A in its real part and imaginary part and does the 
        integral I^y_2 in the notes in the square within the bounds of 
        yminus, yplus, xminus, xplus inputs
    =#
    ans2 = - (dx/8) * (first_back_x(Areal, dx, yplus, xplus)^2 
                    + first_back_x(Aimag, dx, yplus, xplus)^2
                    - first_back_x(Areal, dx, yminus, xplus)^2 
                    - first_back_x(Aimag, dx, yminus, xplus)^2
                    + first_forw_x(Areal, dx, yplus, xminus)^2 
                    + first_forw_x(Aimag, dx, yplus, xminus)^2
                    - first_forw_x(Areal, dx, yminus, xminus)^2 
                    - first_forw_x(Aimag, dx, yminus, xminus)^2 )
    ans2 += - (dx/4) * sum(first_centred_x(Areal, dx, yplus, k)^2 
                    + first_centred_x(Aimag, dx, yplus, k)^2
                    - first_centred_x(Areal, dx, yminus, k)^2 
                    - first_centred_x(Aimag, dx, yminus, k)^2
                    for k in xminus+1:xplus-1)
    return ans2
end


##################################
## Integrals for general thrust ##
##################################

function mod_xth_xflux(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   #=
       Takes a matrix A in its real part and imaginary part and does the 
       integral I^y_2 in the notes in the square within the bounds of 
       yminus, yplus, xminus, xplus inputs

   =#
   ans2 = - (dx/8) * (first_back_x(Areal, dx, yplus, xplus)^2 
                   + first_back_x(Aimag, dx, yplus, xplus)^2
                   - first_forw_x(Areal, dx, yplus, xminus)^2 
                   - first_forw_x(Aimag, dx, yplus, xminus)^2
                   + first_back_x(Areal, dx, yminus, xplus)^2 
                   + first_back_x(Aimag, dx, yminus, xplus)^2
                   - first_forw_x(Areal, dx, yminus, xminus)^2 
                   - first_forw_x(Aimag, dx, yminus, xminus)^2 )
   ans2 += - (dx/4) * sum(first_back_x(Areal, dx, k, xplus)^2 
                   + first_back_x(Aimag, dx, k, xplus)^2
                   - first_forw_x(Areal, dx, k, xminus)^2 
                   - first_forw_x(Aimag, dx, k, xminus)^2
                   for k in yminus+1:yplus-1)
   return ans2
end

function mod_fore_aft(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   ans3 = (1/4)*(-(dx/2)*(Areal[yminus,xplus]^2 + Aimag[yminus,xplus]^2 
                        - Areal[yminus,xminus]^2- Aimag[yminus,xminus]^2 
                        + Areal[yplus,xplus]^2 + Aimag[yplus,xplus]^2 
                        - Areal[yplus,xminus]^2- Aimag[yplus,xminus]^2)
                        -dx * sum(Areal[k,xplus]^2 + Aimag[k,xplus]^2 
                              - Areal[k,xminus]^2 - Aimag[k,xminus]^2 
                                 for k in yminus+1:yplus-1))
   return ans3
end

function mod_thrust_Q(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag)
   # corners
   ans = (1/8)*(dx^2)*(qreal[yplus,xminus]*first_forw_x(hreal,dx,yplus,xminus) +
                        qimag[yplus,xminus]*first_forw_x(himag,dx,yplus,xminus) +
                        qreal[yplus,xplus]*first_back_x(hreal,dx,yplus,xplus) +
                        qimag[yplus,xplus]*first_back_x(himag,dx,yplus,xplus) +
                        qreal[yminus,xminus]*first_forw_x(hreal,dx,yminus,xminus)+
                        qimag[yminus,xminus]*first_forw_x(himag,dx,yminus,xminus)+
                        qreal[yminus,xplus]*first_back_x(hreal,dx,yminus,xplus) +
                        qimag[yminus,xplus]*first_back_x(himag,dx,yminus,xplus) )
   # edges
   ans += (1/4)*(dx^2)*(sum( qreal[k,xminus]*first_forw_x(hreal,dx,k,xminus)
                           + qimag[k,xminus]*first_forw_x(himag,dx,k,xminus)
                           + qreal[k,xplus]*first_back_x(hreal,dx,k,xplus)
                           + qimag[k,xplus]*first_back_x(himag,dx,k,xplus)
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*first_centred_x(hreal,dx,yplus,k)
                           + qimag[yplus,k]*first_centred_x(himag,dx,yplus,k)
                           + qreal[yminus,k]*first_centred_x(hreal,dx,yminus,k)
                           + qimag[yminus,k]*first_centred_x(himag,dx,yminus,k)
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx^2)*(sum( sum( qreal[k,m]*first_centred_x(hreal,dx,k,m)
                                 + qimag[k,m]*first_centred_x(himag,dx,k,m)
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )
   return ans
end 

function mod_cross_term(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   ans4 = (dx/2)*(first_back_x(Areal,dx,yplus,xplus)*first_back_y(Areal,dx,yplus,xplus) 
            + first_back_x(Aimag,dx,yplus,xplus)*first_back_y(Aimag,dx,yplus,xplus)
            - first_back_x(Areal,dx,yminus,xplus)*first_forw_y(Areal,dx,yminus,xplus)
            - first_back_x(Aimag,dx,yminus,xplus)*first_forw_y(Aimag,dx,yminus,xplus))
   ans4 += (dx/2)*(first_forw_x(Areal,dx,yplus,xminus)*first_back_y(Areal,dx,yplus,xminus) 
            + first_forw_x(Aimag,dx,yplus,xminus)*first_back_y(Aimag,dx,yplus,xminus)
            - first_forw_x(Areal,dx,yminus,xminus)*first_forw_y(Areal,dx,yminus,xminus)
            - first_forw_x(Aimag,dx,yminus,xminus)*first_forw_y(Aimag,dx,yminus,xminus))
   ans4 += dx * sum(first_centred_x(Areal,dx,yplus,k)*first_back_y(Areal,dx,yplus,k) 
            + first_centred_x(Aimag,dx,yplus,k)*first_back_y(Aimag,dx,yplus,k)
            - first_centred_x(Areal,dx,yminus,k)*first_forw_y(Areal,dx,yminus,k)
            - first_centred_x(Aimag,dx,yminus,k)*first_forw_y(Aimag,dx,yminus,k)
            for k in xminus+1:xplus-1)

   return (1/2)*ans4
end

function mod_xth_yflux(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   #=
       Takes a matrix A in its real part and imaginary part and does the 
       integral I^y_2 in the notes in the square within the bounds of 
       yminus, yplus, xminus, xplus inputs

   =#
   ans2 = - (dx/8) * (first_back_y(Areal, dx, yplus, xplus)^2 
                   + first_back_y(Aimag, dx, yplus, xplus)^2
                   - first_back_y(Areal, dx, yplus, xminus)^2 
                   - first_back_y(Aimag, dx, yplus, xminus)^2
                   + first_forw_y(Areal, dx, yminus, xplus)^2 
                   + first_forw_y(Aimag, dx, yminus, xplus)^2
                   - first_forw_y(Areal, dx, yminus, xminus)^2 
                   - first_forw_y(Aimag, dx, yminus, xminus)^2 )
   ans2 += - (dx/4) * sum(first_centred_y(Areal, dx, k, xplus)^2 
                   + first_centred_y(Aimag, dx, k, xplus)^2
                   - first_centred_y(Areal, dx, k, xminus)^2 
                   - first_centred_y(Aimag, dx, k, xminus)^2
                   for k in yminus+1:yplus-1)
   return ans2
end


function plane_left_right(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   ans3 = (1/4)*((dx/2)*(Areal[yminus,xplus]^2 + Aimag[yminus,xplus]^2 
                        - Areal[yminus,xminus]^2- Aimag[yminus,xminus]^2 
                        + Areal[yplus,xplus]^2 + Aimag[yplus,xplus]^2 
                        - Areal[yplus,xminus]^2- Aimag[yplus,xminus]^2)
                        +dx * sum(Areal[yplus,k]^2 + Aimag[yplus,k]^2 
                                - Areal[yminus,k]^2 - Aimag[yminus,k]^2 
                                for k in xminus+1:xplus-1))
   return ans3
end

function plane_fore_aft(yminus, yplus, xminus, xplus, dx, Areal, Aimag)
   ans3 = (1/4)*((dx/2)*(Areal[yminus,xplus]^2 + Aimag[yminus,xplus]^2 
                        - Areal[yminus,xminus]^2- Aimag[yminus,xminus]^2 
                        + Areal[yplus,xplus]^2 + Aimag[yplus,xplus]^2 
                        - Areal[yplus,xminus]^2- Aimag[yplus,xminus]^2)
                        + dx * sum(Areal[k,xplus]^2 + Aimag[k,xplus]^2 
                                - Areal[k,xminus]^2 - Aimag[k,xminus]^2 
                                for k in yminus+1:yplus-1))
   return ans3
end

# Adding in a power function
# check the sign on this one
function q_pwr(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag,vel)
   # corners
   ans = (1/8)*(dx^2)*(qreal[yplus,xminus]*himag[yplus,xminus] -
                        qimag[yplus,xminus]*hreal[yplus,xminus] +
                        qreal[yplus,xplus]*himag[yplus,xplus] -
                        qimag[yplus,xplus]*hreal[yplus,xplus] +
                        qreal[yminus,xminus]*himag[yminus,xminus]-
                        qimag[yminus,xminus]*hreal[yminus,xminus]+
                        qreal[yminus,xplus]*himag[yminus,xplus] -
                        qimag[yminus,xplus]*hreal[yminus,xplus] )
   # edges
   ans += (1/4)*(dx^2)*(sum( qreal[k,xminus]*himag[k,xminus]
                           - qimag[k,xminus]*hreal[k,xminus]
                           + qreal[k,xplus]*himag[k,xplus]
                           - qimag[k,xplus]*hreal[k,xplus]
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*himag[yplus,k]
                           - qimag[yplus,k]*hreal[yplus,k]
                           + qreal[yminus,k]*himag[yminus,k]
                           - qimag[yminus,k]*hreal[yminus,k]
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx^2)*(sum( sum( qreal[k,m]*himag[k,m]
                                 - qimag[k,m]*hreal[k,m]
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )

   # add on the velocity term
   ans -= vel*x_thrust_Q(yminus, yplus, xminus, xplus, dx, qreal, qimag,hreal,himag)
   return ans
end 



### Modifications to the expressions to have dx != dy generalisation

# the thrust is a simple generalisation because all the derivatives only depend
# on dx
function sup_x_thrust_Q(yminus, yplus, xminus, xplus, dx, dy, qreal, qimag,hreal,himag)
   # corners
   ans = (1/8)*(dx*dy)*(qreal[yplus,xminus]*first_forw_x(hreal,dx,yplus,xminus) +
                        qimag[yplus,xminus]*first_forw_x(himag,dx,yplus,xminus) +
                        qreal[yplus,xplus]*first_back_x(hreal,dx,yplus,xplus) +
                        qimag[yplus,xplus]*first_back_x(himag,dx,yplus,xplus) +
                        qreal[yminus,xminus]*first_forw_x(hreal,dx,yminus,xminus)+
                        qimag[yminus,xminus]*first_forw_x(himag,dx,yminus,xminus)+
                        qreal[yminus,xplus]*first_back_x(hreal,dx,yminus,xplus) +
                        qimag[yminus,xplus]*first_back_x(himag,dx,yminus,xplus) )
   # edges
   ans += (1/4)*(dx*dy)*(sum( qreal[k,xminus]*first_forw_x(hreal,dx,k,xminus)
                           + qimag[k,xminus]*first_forw_x(himag,dx,k,xminus)
                           + qreal[k,xplus]*first_back_x(hreal,dx,k,xplus)
                           + qimag[k,xplus]*first_back_x(himag,dx,k,xplus)
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*first_centred_x(hreal,dx,yplus,k)
                           + qimag[yplus,k]*first_centred_x(himag,dx,yplus,k)
                           + qreal[yminus,k]*first_centred_x(hreal,dx,yminus,k)
                           + qimag[yminus,k]*first_centred_x(himag,dx,yminus,k)
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx*dy)*(sum( sum( qreal[k,m]*first_centred_x(hreal,dx,k,m)
                                 + qimag[k,m]*first_centred_x(himag,dx,k,m)
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )
   return ans
end 

function sup_q_pwr(yminus, yplus, xminus, xplus, dx, dy, qreal, qimag,hreal,himag,vel)
   # corners
   ans = (1/8)*(dx*dy)*(qreal[yplus,xminus]*himag[yplus,xminus] -
                        qimag[yplus,xminus]*hreal[yplus,xminus] +
                        qreal[yplus,xplus]*himag[yplus,xplus] -
                        qimag[yplus,xplus]*hreal[yplus,xplus] +
                        qreal[yminus,xminus]*himag[yminus,xminus]-
                        qimag[yminus,xminus]*hreal[yminus,xminus]+
                        qreal[yminus,xplus]*himag[yminus,xplus] -
                        qimag[yminus,xplus]*hreal[yminus,xplus] )
   # edges
   ans += (1/4)*(dx*dy)*(sum( qreal[k,xminus]*himag[k,xminus]
                           - qimag[k,xminus]*hreal[k,xminus]
                           + qreal[k,xplus]*himag[k,xplus]
                           - qimag[k,xplus]*hreal[k,xplus]
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*himag[yplus,k]
                           - qimag[yplus,k]*hreal[yplus,k]
                           + qreal[yminus,k]*himag[yminus,k]
                           - qimag[yminus,k]*hreal[yminus,k]
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx*dy)*(sum( sum( qreal[k,m]*himag[k,m]
                                 - qimag[k,m]*hreal[k,m]
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )

   # add on the velocity term
   ans -= vel*sup_x_thrust_Q(yminus, yplus, xminus, xplus, dx, dy, qreal, qimag,hreal,himag)
   return ans
end 

function sup_q_pwr_sign_flip(yminus, yplus, xminus, xplus, dx, dy, qreal, qimag,hreal,himag,vel)
   # corners
   ans = (1/8)*(dx*dy)*(qreal[yplus,xminus]*himag[yplus,xminus] -
                        qimag[yplus,xminus]*hreal[yplus,xminus] +
                        qreal[yplus,xplus]*himag[yplus,xplus] -
                        qimag[yplus,xplus]*hreal[yplus,xplus] +
                        qreal[yminus,xminus]*himag[yminus,xminus]-
                        qimag[yminus,xminus]*hreal[yminus,xminus]+
                        qreal[yminus,xplus]*himag[yminus,xplus] -
                        qimag[yminus,xplus]*hreal[yminus,xplus] )
   # edges
   ans += (1/4)*(dx*dy)*(sum( qreal[k,xminus]*himag[k,xminus]
                           - qimag[k,xminus]*hreal[k,xminus]
                           + qreal[k,xplus]*himag[k,xplus]
                           - qimag[k,xplus]*hreal[k,xplus]
                           for k in yminus+1:yplus-1)
                        +sum(qreal[yplus,k]*himag[yplus,k]
                           - qimag[yplus,k]*hreal[yplus,k]
                           + qreal[yminus,k]*himag[yminus,k]
                           - qimag[yminus,k]*hreal[yminus,k]
                           for k in xminus+1:xplus-1) )
   # bulk
   ans += (1/2)*(dx*dy)*(sum( sum( qreal[k,m]*himag[k,m]
                                 - qimag[k,m]*hreal[k,m]
                                 for m in xminus+1:xplus-1) 
                                 for k in yminus+1:yplus-1) )

   # add on the velocity term
   ans += vel*sup_x_thrust_Q(yminus, yplus, xminus, xplus, dx, dy, qreal, qimag,hreal,himag)
   return ans
end 

