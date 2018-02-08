
#int nb_sim = number of simulations
#matrix float dim x 2 boundaries = matrix of boundaries, size dim x 2 (for each dim, provides min and max)
#float time = end time of simulation

function simulations_CPU(nb_sim, boundaries, time, model, sig_model, dim)

    end_state=zeros(nb_sim,dim)
    init_state=zeros(dim)

    for i = 1:nb_sim

        for j=1:dim
            init_state[j]=rand([boundaries[1,j],boundaries[2,j]])
        end

        prob_sde = SDEProblem(model,sig_model,init_state,(0.0,time))
        sol = solve(prob_sde)

        tfinal=size(sol.u)[1]

        for k=1:dim
            end_state[i,k]=sol.u[tfinal][k]
        end
    end

    return end_state
end
