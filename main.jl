using Distributions, Plots, DifferentialEquations


##### DEFINITION OF THE PROBLEM
nb_bin=20 #discretization of the space in each dimension for the frequency matrix
nb_sim=1000 #number of simmulations
end_time= 10.0 #end time for each simulation


##### DEFINITION OF THE MODEL
a1=1
a2=1
b1=1
b2=1
k1=1
k2=1
n=4
S=0.5
lambda=0.01

##### 2D model
dim=2
boundaries=[[0.0,3.0] [0.0,3.0]] #boundaries for each dimension
tuple_dim=(nb_bin,nb_bin) #tuple of dimension of the discretized space
function model(du,u,p,t)
    du[1] = a1*u[1]^n/(u[1]^n+S^n)+b1*S^n/(u[2]^n+S^n)-k1*u[1]
    du[2] = a2*u[2]^n/(u[2]^n+S^n)+b2*S^n/(u[1]^n+S^n)-k2*u[2]
end

function sig_model(du,u,p,t)
    du[1] = sqrt(abs(a1*u[1]^n/(u[1]^n+S^n)+b1*S^n/(u[2]^n+S^n)-k1*u[1]))
    du[2] = sqrt(abs(a2*u[2]^n/(u[2]^n+S^n)+b2*S^n/(u[1]^n+S^n)-k2*u[2]))
end

##### 3D model
dim=3 #dimension of the network
boundaries=[[0.0, 1.5] [0,0.4] [1.4, 1.5]] #boundaries for each dimension
tuple_dim=(nb_bin,nb_bin,nb_bin) #tuple of dimension of the discretized space
function model(du,u,p,t)
    du[1] = u[3]*u[1]^n/(u[1]^n+S^n)+b1*S^n/(u[2]^n+S^n)-k1*u[1]
    du[2] = u[3]*u[2]^n/(u[2]^n+S^n)+b2*S^n/(u[1]^n+S^n)-k2*u[2]
    du[3] = -lambda*u[3]
end

function sig_model(du,u,p,t)
    du[1] = sqrt(abs(u[3]*u[1]^n/(u[1]^n+S^n)+b1*S^n/(u[2]^n+S^n)-k1*u[1]))
    du[2] = sqrt(abs(u[3]*u[2]^n/(u[2]^n+S^n)+b2*S^n/(u[1]^n+S^n)-k2*u[2]))
    du[3] = sqrt(abs(-lambda*u[3]))
end


##### SIMULATIONS
state=simulations_CPU(nb_sim,boundaries,end_time,model,sig_model,dim)

##### PFM
U=Prob_Flux_Method(state,boundaries,tuple_dim)
span=Span(boundaries,nb_bin,dim)


##### PFM reduction 2D
U_pfm=PFM_reduction_2D(state,boundaries,nb_bin,1,2)
span_2D=Span_2D(boundaries,nb_bin,1,2)


##### KDE reduction 2D
U_kde=KDE_reduction_2D(state,boundaries,nb_bin,1,2)
span_2D=Span_2D(boundaries,nb_bin,1,2)

x=span_2D[1,:]
y=span_2D[2,:]

surface(y,x,U_pfm)
