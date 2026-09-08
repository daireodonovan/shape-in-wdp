# Storage for the finite difference functions that I use in other files

function first_centred_y(A,dx,y,x)
    #=
      Pass a matrix as an argument and then calculate the first order derivative
      using a centred difference approximation given the step size (dx) and the 
      the point in the array (y,x).
    =#  
    return (-A[y-1,x] + A[y+1,x])/(2*dx)
end

function first_back_y(A,dx,y,x)
    #=
      Pass a matrix as an argument and then calculate the first order derivative
      using a backward difference given the step size (dx) and the 
      the point in the array (y,x).
    =#  
    return (3*A[y,x] - 4*A[y-1,x] + A[y-2,x])/(2*dx)
end

function first_forw_y(A,dx,y,x)
    #=
      Pass a matrix as an argument and then calculate the first order derivative
      using a forward difference given the step size (dx) and the 
      the point in the array (y,x).
    =#  
    return (-3*A[y,x] + 4*A[y+1,x] - A[y+2,x])/(2*dx)
end

function first_centred_x(A,dx,y,x)
  #=
    Pass a matrix as an argument and then calculate the first order derivative
    using a centred difference approximation given the step size (dx) and the 
    the point in the array (y,x).
  =#  
  return (-A[y,x-1] + A[y,x+1])/(2*dx)
end

function first_back_x(A,dx,y,x)
  #=
    Pass a matrix as an argument and then calculate the first order derivative
    using a backward difference given the step size (dx) and the 
    the point in the array (y,x).
  =#  
  return (3*A[y,x] - 4*A[y,x-1] + A[y,x-2])/(2*dx)
end

"""
    Pass a matrix as an argument and then calculate the first order derivative
    using a forward difference given the step size (dx) and the 
    the point in the array (y,x).
"""  
function first_forw_x(A,dx,y,x)

  return (-3*A[y,x] + 4*A[y,x+1] - A[y,x+2])/(2*dx)
end

##################
## Second order ##
##################

# x 

function second_centred_x(A,dx,y,x)
  return (A[y,x-1] - 2*A[y,x] + A[y,x+1])/(dx^2)  
end

function second_forw_x(A,dx,y,x)
  return (2*A[y,x] - 5*A[y,x+1] + 4*A[y,x+2] -A[y,x+3]) / (dx^2)
end

function second_back_x(A,dx,y,x)
  return (2*A[y,x] - 5*A[y,x-1] + 4*A[y,x-2] -A[y,x-3]) / (dx^2)
end

# y

function second_centred_y(A,dx,y,x)
  return (A[y-1,x] - 2*A[y,x] + A[y+1,x])/(dx^2)  
end

function second_forw_y(A,dx,y,x)
  return (2*A[y,x] - 5*A[y+1,x] + 4*A[y+2,x] -A[y+3,x]) / (dx^2)
end

function second_back_y(A,dx,y,x)
  return (2*A[y,x] - 5*A[y-1,x] + 4*A[y-2,x] -A[y-3,x]) / (dx^2)
end

#= 
  Higher order approximations
=# 

# order 4
function first_forw_x_ho4(A,dx,y,x)

  return ((-25*A[y,x] + 48*A[y,x+1] - 36*A[y,x+2]+ 16*A[y,x+3] - 3*A[y,x+4])/(12*dx))
end

function first_back_x_ho4(A,dx,y,x)

  return ((25*A[y,x] - 48*A[y,x-1] + 36*A[y,x-2]- 16*A[y,x-3] + 3*A[y,x-4])/(12*dx))
end

function first_centred_x_ho4(A,dx,y,x)

  return ((A[y,x-2] - 8*A[y,x-1] + 8*A[y,x+1] - A[y,x+2])/(12*dx))
end



function first_forw_y_ho4(A,dx,y,x)

  return ((-25*A[y,x] + 48*A[y+1,x] - 36*A[y+2,x]+ 16*A[y+3,x] - 3*A[y+4,x])/(12*dx))
end

function first_back_y_ho4(A,dx,y,x)

  return ((25*A[y,x] - 48*A[y-1,x] + 36*A[y-2,x]- 16*A[y-3,x] + 3*A[y-4,x])/(12*dx))
end

function first_centred_y_ho4(A,dx,y,x)

  return ((A[y-2,x] - 8*A[y-1,x] + 8*A[y+1,x] - A[y+2,x] )/(12*dx))
end


function second_centred_x_ho4(A,dx,y,x)
  return ((-A[y,x-2] + 16*A[y,x-1] -30*A[y,x] +16*A[y,x+1] -A[y,x+2])/(12*dx^2) )
end

function second_forw_x_ho4(A,dx,y,x)
  return ((45*A[y,x] - 154*A[y,x+1] + 214*A[y,x+2] - 156*A[y,x+3] +61*A[y,x+4] 
        - 10*A[y,x+5])/(12*dx^2)  )
end

function second_back_x_ho4(A,dx,y,x)
  return ((45*A[y,x] - 154*A[y,x-1] + 214*A[y,x-2] - 156*A[y,x-3] +61*A[y,x-4] 
        - 10*A[y,x-5] )/(12*dx^2) )
end


function second_centred_y_ho4(A,dx,y,x)
  return ((-A[y-2,x] + 16*A[y-1,x] - 30*A[y,x] + 16*A[y+1,x]-A[y+2,x])/(12*dx^2)  )
end

function second_forw_y_ho4(A,dx,y,x)
  return ((45*A[y,x] - 154*A[y+1,x] + 214*A[y+2,x] - 156*A[y+3,x] +61*A[y+4,x] 
        - 10*A[y+5,x])/(12*dx^2)  )
end

function second_back_y_ho4(A,dx,y,x)
  return ((45*A[y,x] - 154*A[y-1,x] + 214*A[y-2,x] - 156*A[y-3,x] +61*A[y-4,x] 
        - 10*A[y-5,x])/(12*dx^2)  )
end

# order 6

function first_forw_x_ho6(A,dx,y,x)

  return ((-147*A[y,x] + 360*A[y,x+1] - 450*A[y,x+2]+ 400*A[y,x+3] - 225*A[y,x+4]
                        +72*A[y,x+5] - 10*A[y,x+6])/(60*dx))
end

function first_back_x_ho6(A,dx,y,x)

  return ((147*A[y,x] - 360*A[y,x-1] + 450*A[y,x-2]- 400*A[y,x-3] + 225*A[y,x-4]
                        -72*A[y,x-5] + 10*A[y,x-6])/(60*dx))
end

function first_centred_x_ho6(A,dx,y,x)

  return ((-A[y,x-3] + 9*A[y,x-2] - 45*A[y,x-1] + 45*A[y,x+1] - 9*A[y,x+2] 
          + A[y,x+3]  )/(60*dx))
end



function first_forw_y_ho6(A,dx,y,x)

  return ((-147*A[y,x] + 360*A[y+1,x] - 450*A[y+2,x]+ 400*A[y+3,x] - 225*A[y+4,x]
                        +72*A[y+5,x] - 10*A[y+6,x])/(60*dx))
end

function first_back_y_ho6(A,dx,y,x)

  return ((147*A[y,x] - 360*A[y-1,x] + 450*A[y-2,x]- 400*A[y-3,x] + 225*A[y-4,x]
                        -72*A[y-5,x] + 10*A[y-6,x])/(60*dx))
end

function first_centred_y_ho6(A,dx,y,x)

  return ((-A[y-3,x] + 9*A[y-2,x] - 45*A[y-1,x] + 45*A[y+1,x] - 9*A[y+2,x] 
          + A[y+3,x]  )/(60*dx))
end


function second_centred_x_ho6(A,dx,y,x)
  return ((A[y,x-3] - 13.5*A[y,x-2] + 135*A[y,x-1] - 245*A[y,x] 
          + 135*A[y,x+1]-13.5*A[y,x+2]+A[y,x+3])/(90*dx^2)  )
end

function second_forw_x_ho6(A,dx,y,x)
  return ((938*A[y,x] - 4014*A[y,x+1] + 7911*A[y,x+2] - 9490*A[y,x+3] +7380*A[y,x+4] 
        - 3618*A[y,x+5] + 1019*A[y,x+6] - 126*A[y,x+7])/(180*dx^2)  )
end

function second_back_x_ho6(A,dx,y,x)
  return ((938*A[y,x] - 4014*A[y,x-1] + 7911*A[y,x-2] - 9490*A[y,x-3] +7380*A[y,x-4] 
        - 3618*A[y,x-5] + 1019*A[y,x-6] - 126*A[y,x-7])/(180*dx^2)  )
end


function second_centred_y_ho6(A,dx,y,x)
  return ((A[y-3,x] - 13.5*A[y-2,x] + 135*A[y-1,x] - 245*A[y,x] 
          + 135*A[y+1,x]-13.5*A[y+2,x]+A[y+3,x])/(90*dx^2)  )
end

function second_forw_y_ho6(A,dx,y,x)
  return ((938*A[y,x] - 4014*A[y+1,x] + 7911*A[y+2,x] - 9490*A[y+3,x] +7380*A[y+4,x] 
        - 3618*A[y+5,x] + 1019*A[y+6,x] - 126*A[y+7,x])/(180*dx^2)  )
end

function second_back_y_ho6(A,dx,y,x)
  return ((938*A[y,x] - 4014*A[y-1,x] + 7911*A[y-2,x] - 9490*A[y-3,x] +7380*A[y-4,x] 
        - 3618*A[y-5,x] + 1019*A[y-6,x] - 126*A[y-7,x])/(180*dx^2)  )
end