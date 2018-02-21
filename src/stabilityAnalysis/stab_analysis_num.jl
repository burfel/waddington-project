
###########################################################################
###################### Numerial stability analysis ########################
###########################################################################

using SymEngine, SymPy, DifferentialEquations

# For range = 1:5
range_inf = 0;
range_sup = 5;


function stab_analysis(model, np_model, u0, tspan, p, range_inf, range_sup)


    function string_as_varname(s::String,v::Any)
            s = Symbol(s)
            @eval (($s) = ($v))
    end

    # Convert expressions of the parameters and variables into SymPy elements in order to later evaluate mode.funcs and so that they are converted into SymPy expressions.

    for i in 1:length(model.params)
        string_as_varname(string(model.params[i]), p[i])
    end



        ###### Stop considering them as parameters in the parameterized model #########



        # Define ST problem and solve it (note: parameters are then automatically substituted!):

        fixed_points = [];
        numb_fp = 0;

        global fixed_points

        for i in 1:100   # We get an arrays of arrays (one per sol. vector).

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

        # Convert model.params and model.syms (Expr) to SymEngine type

        for counter in 1:length(model.params)

            string_as_varname(string(model.params[counter]), SymEngine.symbols(model.params[counter]))
        end

        for counter in 1:length(model.syms)

            string_as_varname(string(model.syms[counter]), SymEngine.symbols(model.syms[counter]))
        end

        # Create dictionary of parameters and array of dictionaries of the variables (one dic per fixed point).

        par_dic = Dict(zip([eval.(model.params)[i] for i in 1:length(model.params)], [p[i] for i in 1:length(model.params)]));
        var_dic_set = [];

        for counter in 1:length(fixed_points) # (Creates arrays of dictionaries of the form SymEngine => Float64).

                push!(var_dic_set, Dict(zip([eval(model.syms[i]) for i in 1:length(model.syms)], [ (fixed_points[counter])[i] for i in 1:length(model.syms)])))
        end

        # Substitute Jacobian with parameters

        J_par = J;

        for counter in 1:length(model.params)

            J_par = SymEngine.subs.(J_par, eval(model.params[counter]) => par_dic[eval(model.params[counter])] )
        end


        # Substitute jacobian [SymEngine type] with fixed points and create matrix of jacobians (one for each fixed point).

        JM = [];

        Len = length(fixed_points[1]); # Necessary because if only one sol. the result is not an array but a column vector.


        for counter1 in 1:numb_fp

            J_par_local = J_par

            for counter2 in 1:length(model.syms)

                J_par_local = SymEngine.subs.(J_par_local, eval(model.syms[counter2]) => ( (var_dic_set[counter1])[eval(model.syms[counter2])] ) )
            end

            push!(JM, J_par_local)
        end

        # Convert JM matrices to numeric matrices.

        JM = Array{Array{Float64, 2}}(JM);



        ######### Stability Analysis ##########



        # Evaluate eigenvalues of each Jacobian and store fixed points in vectors (whose elements, that is, each fixed point, is in the form of a dictionary):

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

        println("stable = ", round.(stable,1))
        println("unstable = ", round.(unstable,1))
        println("saddle = ", round.(saddle,1))
        println("imag = ", round.(imag,1))


    end
