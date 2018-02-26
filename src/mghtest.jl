# -*- coding: utf-8 -*-
#=

Function that plots the tries out different plotting methods to visualise ODE/ SDE solutions, contour maps.

Note: -

Args:
    -

Returns:
    -

=#

using Plots
plotly()
plot(rand(5,5),linewidth=2,title="Epigenetic landscape")

using DifferentialEquations
f(u,p,t) = 1.01*u
u0=1/2
tspan = (0.0,1.0)

# SOLUTION OF ODE PROBLEM
prob_ode = ODEProblem(f,u0,tspan)
sol_ode = solve(prob_ode)
# PLOT OF SOLUTION
plot(title = "Solution of the ODE problem: ", sol_ode)

function lorenz(du,u,p,t)
  du[1] = 10.0(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end

function σ_lorenz(du,u,p,t)
  du[1] = 3.0
  du[2] = 3.0
  du[3] = 3.0
end

# SOLUTION OF SDE PROBLEM WITH DIAGONAL NOISE
prob_sde_lorenz = SDEProblem(lorenz,σ_lorenz,[1.0,0.0,0.0],(0.0,10.0))
sol_sde = solve(prob_sde_lorenz)
# PLOT OF SOLUTION
plot(title = "Solution of the SDE problem: ", sol_sde,vars=(1,2,3))

##############################################################################

# SOME RANDOM CONTOUR MAPS
x = 1:0.5:20
y = 1:0.5:10
f(x,y) = begin
        (3x + y ^ 2) * abs(sin(x) + cos(y))
    end
X = repmat(x',length(y),1)
Y = repmat(y,1,length(x))
Z = map(f,X,Y)
p1 = contour(x,y,f,fill=true)
p2 = contour(x,y,Z)
plot(p1,p2)

# TODO:  would be nice to add this plot: https://plot.ly/julia/contour-plots/
# also: plots from Gadfly.jl, which isn't working currently.

#=
using MultivariateStats
X = rand(50,50);
M = fit(PPCA, X, maxoutdim=2)
Y = transform(M,X)
XR = reconstruct(M,Y)

X = rand(50,50);
M = fit(KernelPCA,X;maxoutdim=10, inverse = true)

using DataStructures
using ManifoldLearning
X = rand(50,50)
Y = transform(LEM, X; k=12, d = 2, t = 1.0)
=#
