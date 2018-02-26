
###########################################################################
###################### Numerial stability analysis ########################
###########################################################################
# -*- coding: utf-8 -*-
#=

Program that performs stability analysis.

Note: -

Args:
    model:
    np_model:
    u0:
    tspan:
    p:
    range_inf (optional):
    range_sup (optional):

Returns:
    -

=#

using SymEngine, SymPy, DifferentialEquations

function stab_analysis(model, np_model, u0, tspan, p, range_inf=0, range_sup=5)

    # Define function that converts a variable into a symbol, so
    # you can assign different values to it at different times;
    # handy to evaluate model functions at different values.
    function string_as_varname(s::String,v::Any)
            s = Symbol(s)
            @eval (($s) = ($v))
    end

    # Convert model parameters into SymPy elements
    # in order to later evaluate mode.funcs
    for i in 1:length(model.params)
        string_as_varname(string(model.params[i]), p[i])
    end

    # Define SDE problem and solve it
    # Note: parameters are then automatically substituted

    fixed_points = [];
    numb_fp = 0;

    global fixed_points

    for i in 1:100   # ...we get an array of arrays (one per solution vector).

        u0 = rand(length(model.syms))*(range_sup-range_inf) - range_inf
        prob = SteadyStateProblem(np_model,u0)
        sol = DifferentialEquations.solve(prob,SSRootfind())

        if !in(round.(sol,1),round.(fixed_points, 1)) && all(x->x >= 0, sol)
            push!(fixed_points, sol)
            numb_fp += 1;
        end
    end

    ##### Now change type of all the free variables from SymPy to SymEngine ######

    # Jacobian [SymEngine type]

    J = model.symjac;

    # Converts model.params and model.syms (Expr) to SymEngine type
    for counter in 1:length(model.params)
        string_as_varname(string(model.params[counter]), SymEngine.symbols(model.params[counter]))
    end

    for counter in 1:length(model.syms)
        string_as_varname(string(model.syms[counter]), SymEngine.symbols(model.syms[counter]))
    end

    # Creates a dictionary of parameters and
    # and an array of dictionaries of the variables (one dic per fixed point).
    par_dic = Dict(zip([eval.(model.params)[i] for i in 1:length(model.params)], [p[i] for i in 1:length(model.params)]));
    var_dic_set = [];

    # Creates arrays of dictionaries of the form SymEngine => Float64
    for counter in 1:length(fixed_points)
            push!(var_dic_set, Dict(zip([eval(model.syms[i]) for i in 1:length(model.syms)], [ (fixed_points[counter])[i] for i in 1:length(model.syms)])))
    end

    # Substitutes Jacobian with parameters
    J_par = J;

    for counter in 1:length(model.params)
        J_par = SymEngine.subs.(J_par, eval(model.params[counter]) => par_dic[eval(model.params[counter])] )
    end


    # Substitutes Jacobian [SymEngine type] with fixed points and create matrix of jacobians (one for each fixed point).
    JM = [];
    Len = length(fixed_points[1]); # Necessary because if only one sol. the result is not an array but a column vector.

    for counter1 in 1:numb_fp
        J_par_local = J_par

        for counter2 in 1:length(model.syms)
            J_par_local = SymEngine.subs.(J_par_local, eval(model.syms[counter2]) => ( (var_dic_set[counter1])[eval(model.syms[counter2])] ) )
        end

        push!(JM, J_par_local)
    end

    # Converts JM matrices to numeric matrices.
    JM = Array{Array{Float64, 2}}(JM);


    ######### Actual stability analysis ##########

    # Evaluate eigenvalues of each Jacobian and store fixed points in vectors
    # whose elements, ie each fixed point, is in the form of a dictionary)

    stable = [];
    saddle = [];
    imag = [];
    unstable = [];

    eig_vals = [];
    global eig_vals

    for i in 1:length(JM)

        eig = eigvals(JM[i]);
        push!(eig_vals, eig)

        if (real(eig).<0) == trues(similar(eig)) # If all eig < 0
            push!(stable, fixed_points[i])

        elseif (real(eig).>0) == trues(similar(eig)) # If all eig > 0
            push!(unstable, fixed_points[i])

        elseif (real(eig).==0) == trues(similar(eig)) # If only imaginary terms
            push!(imag, fixed_points[i])

        else
            push!(saddle, fixed_points[i])

        end
    end

    println("In range ", range_inf, " to ", range_sup, ":")

    println("stable: ", round.(stable,1))
    println("unstable: ", round.(unstable,1))
    println("saddle: ", round.(saddle,1))
    println("imag: ", round.(imag,1))


end

#-----TODO: Call function/ Provide sample code.
#stab_analysis(model, np_model, u0, tspan, p, range_inf=0, range_sup=5)
#
