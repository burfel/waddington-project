using Plots
plotly()
plot(rand(5,5),linewidth=2,title="My Plot")#


using DifferentialEquations
f(u,p,t) = 1.01*u
u0=1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
plot(sol)

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

prob_sde_lorenz = SDEProblem(lorenz,σ_lorenz,[1.0,0.0,0.0],(0.0,10.0))
sol = solve(prob_sde_lorenz)
plot(sol,vars=(1,2,3))

##############################################################################


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


# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a kernel PCA model
#M = fit(KernelPCA, Xtr; maxoutdim=100, inverse=true)
