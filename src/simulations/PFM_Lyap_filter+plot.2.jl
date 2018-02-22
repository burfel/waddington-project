
# Definition of the required functions:

function sq_region(min=0.0, max=3.0, span=20)

    dimx1 = collect(linspace(min, max, span))
    dimx2 = dimx1
    coordinates_matrix = []; # Coordinates matrix
    for n in dimx2
           for i in dimx1
               push!(coordinates_matrix, [i,n])
           end
    end

    coordinates_matrix = reshape(coordinates_matrix, length(dimx1),:)
    global coordinates_matrix = coordinates_matrix[1:(end-1),1:(end-1)]
end

function approx_grad_2Dvec(M) # Calculate approximate gradient and adjust dimensions if M == array of 2D vectors.

         dMx = (diff(M,1))[:,1:(size(diff(M,1))[1])]
         dMy = (diff(M,2))[1:(size(diff(M,2))[2]),:]

         grad_M = [[dMx[i][1],dMy[i][2]] for i in 1:length(dMx)]
         grad_M = reshape(grad_M, length(dMx[:,1]), :) # 19x19

         return grad_M
end

function approx_grad_scalar_in_2D(M) # Calculate approximate gradient and adjust dimensions if M == array of values. Output: array of 2D vectors.

         dMx = (diff(M,1))[:,1:(size(diff(M,1))[1])]
         dMy = (diff(M,2))[1:(size(diff(M,2))[2]),:]

         grad_M = [[dMx[i],dMy[i]] for i in 1:length(dMx)]
         grad_M = reshape(grad_M, length(dMx[:,1]), :) # 19x19

         return grad_M
end

function f_grad(model, p, coordinates_matrix)   # This function gets the function from a parameterized ODE system, first substitutes the
                                                # parameters, and finally each variable with the corresponding coordinates.
                                                # Note: coordinates_matrix should be an array containing the coordinates of every U point.
    function string_as_varname(s::String,v::Any)
            s = Symbol(s)
            @eval (($s) = ($v))
    end

    # Convert expressions of the parameters and variables into SymPy elements in order to later evaluate mode.funcs and so that they are converted into SymPy expressions.

    for i in 1:length(model.params)
        string_as_varname(string(model.params[i]), p[i])
    end

    # Get the value of f at each point (remember that f is a vector): this creates an array of f vectors, F.

    global F = [];

    for c in coordinates_matrix
        for i in 1:length(model.syms)
            string_as_varname(string(model.syms[i]), c[i])
        end
        push!(F, eval.(model.funcs))
    end

    F = (reshape(F, length(coordinates_matrix[:,1]),:)) # 19x19

end

function filter(v_field)

    global non_Lyap = [];
    global Lyap = [];

    L = length(v_field[1])

    for i in 1:length(v_field)
        if (v_field[i])[3] < 0
            push!(Lyap, [i, v_field[i]])
        else
            push!(non_Lyap, [i, v_field[i]])
        end
    end
end

function filter_scatter(v_field)

    global non_Lyap2 = [];
    global Lyap2 = [];

    L = length(v_field[1])

    for i in 1:length(v_field)
        if (v_field[i])[3] < 0
            push!(Lyap2, U_18[i])
        else
            push!(non_Lyap2, U_18[i])
        end
    end
end


function fplot(U_19, Lyap, non_Lyap, Ly, non_Ly)

    U_18 = U_19[1:end-1, 1:end-1]

    if Ly==1

     # Since this is the z-axis coordinates, the size needs to fit the derivative-matrices sizes.

    global x_coord_Lyap = [coordinates_matrix[i[1]][1] for i in Lyap]
    global y_coord_Lyap = [coordinates_matrix[i[1]][2] for i in Lyap]

    z_coord_Lyap = [];

    for i in Lyap
        push!(z_coord_Lyap, U_18[i[1]])
    end

    Lyap_x = [i[2][1] for i in Lyap]
    Lyap_y = [i[2][2] for i in Lyap]
    Lyap_z = [i[2][3] for i in Lyap]

    pygui(true)
    fig = figure()
    ax = fig[:gca](projection="3d")

    ax[:quiver](x_coord_Lyap, y_coord_Lyap, z_coord_Lyap, Lyap_x, Lyap_y, Lyap_z, color="b", LineWidth=0.3)
end
if non_Ly==1

    # Now plot the non_Lyap vectors on the top of the Lyap ones:

    global x_coord_non_Lyap = [coordinates_matrix[i[1]][1] for i in non_Lyap]
    global y_coord_non_Lyap = [coordinates_matrix[i[1]][2] for i in non_Lyap]

    z_coord_non_Lyap = [];

    for i in non_Lyap
        push!(z_coord_non_Lyap, U_18[i[1]])
    end

    non_Lyap_x = [i[2][1] for i in non_Lyap]
    non_Lyap_y = [i[2][2] for i in non_Lyap]
    non_Lyap_z = [i[2][3] for i in non_Lyap]

    pygui(true)
    fig = figure()
    ax = fig[:gca](projection="3d")
    ax[:quiver](x_coord_non_Lyap, y_coord_non_Lyap, z_coord_non_Lyap, non_Lyap_x, non_Lyap_y, non_Lyap_z, color="g", LineWidth=0.3)
end
end


function vfield(F, grad_U, dFc_vec; len=1)

    grad_U_18 = grad_U[1:end-1, 1:end-1] # 18x18 (now grad_U and F need to be adjusted to fit the size of dFc)

    dU_dx = reshape([[grad_U_18[i][1], 0] for i in 1:length(grad_U_18)], size(grad_U_18)[1], :)
    dU_dy = dU_dx = reshape([[0, grad_U_18[i][1]] for i in 1:length(grad_U_18)], size(grad_U_18)[1], :)
    D = abs.(F)[1:end-1, 1:end-1]

    dU_dt = reshape([-dU_dx[i].^2 - dU_dy[i].^2 + vecdot(dFc_vec[i],D[i]) for i in 1:length(dU_dx)], size(dU_dx)[1], :) # 18x18

    dx_dt = [i[1] for i in F][1:end-1, 1:end-1]
    dy_dt = [i[2] for i in F][1:end-1, 1:end-1]


    # Vectorize by setting (dx_dt, dy_dt, dU_dt)

    global v_field = normalize.(reshape([[dx_dt[i], dy_dt[i], (dU_dt[i][1] + dU_dt[i][2])] for i in 1:length(dU_dt)], size(dU_dt)[1], :))./len

end

###############################################################################################
####### Get vector field (dx/dt, dy/dt, dU/dt)*(dt, dt, dt), filter it as Lyap (dU/dt <0) #####
#######  or non_Lyap (>0), and plot those 2 classes in different colours (surface plot)   #####
#######  // plot (x, y, U), also differentiating Lyap and non_Lyap (scatter plot)         #####
###############################################################################################

# Enter the potential matrix and a matrix containing their coordinates:

# span = 20
# U = rand(span,span) # Potential matrix

sq_region(min=0.0, max=3.0, span=20)

# Now, curl = grad_U + F, so first calculate grad_U:

grad_U = approx_grad_scalar_in_2D(U)

# Secondly, calculate F by subsituting the model function in each coordinate. Get the function and substitute with symbolic algebra.

## (Example input:)

# model = @ode_def toymodel begin
#     dx1 = a*x1^n/(x1^n+S^n)+b1*S^n/(x2^n+S^n)-k1*x1
#     dx2 = a*x2^n/(x2^n+S^n)+b2*S^n/(x1^n+S^n)-k2*x2
# end a b1 b2 k1 k2 n S

# u0=zeros(2); tspan=(0.0,10.0); p=[1,1,1,1,1,4,0.5];

f_grad(model, p, coordinates_matrix)


# Finally, get the 'curl' force

Fc = grad_U + F

# Take the approximate gradient of the 'curl' force

dFc_vec = approx_grad_2Dvec(Fc) # Fc is 19x19, so we get m-1 x n-1.


# Calculate dx/dt, dy/dt, and dU/dt

vfield(F, grad_U, dFc_vec)

# if wanted to try with differentials instead of derivatives:
#dt = 0.3
# v_field_diff = reshape([[dx_dt[i]*dt, dy_dt[i]*dt, (dU_dt[i][1] + dU_dt[i][2])*dt] for i in 1:length(dU_dt)], size(dU_dt)[1], :) # This is (dx, dy, dU).

# Filter the f field in 3D by wether they point up while the grad_U in 3D points down.

filter(v_field)

filter_scatter(v_field)

scatter(x_coord_non_Lyap, y_coord_non_Lyap, non_Lyap2)
scatter(x_coord_Lyap, y_coord_Lyap, Lyap2)



# Plot non_Lyap part as default with quiver:

fplot(U_19, Lyap, non_Lyap, Ly, non_Ly)
