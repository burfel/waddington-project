using Distributions
using Plots
using DifferentialEquations


#------DEFINITION OF THE PROBLEM------

# dimension of the network
dim = 2 

# boundaries for each dimension
boundaries = [[0., 3.] [0. ,3.]] 

# discretisation of the space in each dimension for the frequency matrix
nb_bin = 20 

# tuple of dimension of the discretised space
tuple_dim = (nb_bin, nb_bin) 

# number of simulations
nb_sim = 10000 

# end time for each simulation
end_time = 10. 


#------DEFINITION OF THE MODEL--------
a1 = 1
a2 = 1
b1 = 1
b2 = 1
k1 = 1
k2 = 1
n = 4
S = .5
lambda = .01

function model(du, u, p ,t)
    du[1] = a1 * u[1]^n / (u[1]^n + S^n) + b1 * S^n / (u[2]^n+S^n) - k1 * u[1]
    du[2] = a2 * u[2]^n / (u[2]^n + S^n) + b2 * S^n / (u[1]^n+S^n) - k2 * u[2]
end

function sig_model(du, u, p, t)
    du[1] = sqrt(abs(a1 * u[1]^n / (u[1]^n + S^n) + b1 * S^n / (u[2]^n+S^n) - k1 * u[1]))
    du[2] = sqrt(abs(a2 * u[2]^n / (u[2]^n + S^n) + b2 * S^n / (u[1]^n+S^n) - k2 * u[2]))
end


#------SIMULATIONS---------
state = simulations_CPU(nb_sim, boundaries, end_time, model, sig_model, dim)


#------Probability flux method (PFM)--------
U = Prob_Flux_Method(state, boundaries, tuple_dim)
span = Span(boundaries, nb_bin, dim)


#------PFM reduction 2D------------
U_2D = PFM_reduction_2D(state, boundaries, nb_bin, 1, 2)
span_2D = Span_2D(boundaries, nb_bin, 1, 2)

x = span_2D[1,:]
y = span_2D[2,:]

surface(x, y, U_2D)
