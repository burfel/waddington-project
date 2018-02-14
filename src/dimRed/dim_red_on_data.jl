using MultivariateStats
using DataStructures
data_array = convert(Array, data)
#data_array = collect(skipmissing(data_array))
data_array = reshape(data_array, (547,96))



M1 = fit(PCA, data_array; maxoutdim=2)
Y1 = transform(M1, data_array)
X_PCA = reconstruct(M1, Y1)


M2 = fit(PPCA, data_array; maxoutdim=2)
Y2 = transform(M2, data_array)
X_PPCA = reconstruct(M2, Y2)


using Plots

scatter(data_array)


scatter(Y1,title="Y1 - PCA transform")
scatter(X_PCA,title="X_PCA - PCA reconstruct")

scatter(Y2,title="Y2 - PPCA transform")
scatter(X_PPCA,title="X_PPCA - PPCA reconstruct")


using ManifoldLearning

# Isomap
Y_Isomap = transform(Isomap,data_array; k=12,d=2)

# Diffusion maps
Y_DiffMap = transform(DiffMap, data_array; d=2, t=1, ɛ=1.0)

# Laplacian Eigenmaps -- working
Y_LEM = transform(LEM,data_array; k=12, d=2, t=1.0)

# Local linear embedding
Y_LLE = transform(LLE, data_array; k = 12, d = 2)

# Hessian Eigenmaps transformation
Y_HLLE = transform(HLLE, data_array; k = 12, d = 2)

# Local tangent space alignment -- working
Y_LTSA = transform(LTSA, data_array; k = 12, d = 2)

scatter(Y_LEM, title="LEM - ManifoldLearning")

scatter(Y_LTSA, title="LTSA - ManifoldLearning")
