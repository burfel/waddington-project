using MultivariateStats
using DataStructures
data_array = convert(Array,data)
data_array = collect(skipmissing(data_array))
data_array = reshape(data_array,(547,96))
data_array = data_array'


M1 = fit(PCA, data_array;maxoutdim=2)
Y1 = transform(M1,data_array)
X_PCA = reconstruct(M1,Y1)
projection(X_PCA)

M2 = fit(PPCA, data_array; maxoutdim=2)
Y2 = transform(M2,data_array)
X_PPCA = reconstruct(M2,Y2)

# train a kernel PCA model
M3 = fit(KernelPCA, data_array; maxoutdim=2, inverse=true)
# apply kernel PCA model to testing set
Y3 = transform(M3, data_array)
# reconstruct testing observations (approximately)
X_KPCA = reconstruct(M3, data_array)

M4 = fit(ICA, data_array,k=2)
Y4 = transform(M4,data_array)


# train a FactorAnalysis model
M5 = fit(FactorAnalysis, data_array; maxoutdim=2)
# apply FactorAnalysis model to testing set
Y5 = transform(M5, data_array)
# reconstruct testing observations (approximately)
X_FA = reconstruct(M5, Y5)

using Plots

pdata = scatter(data_array);


scatter(Y1,title="Y1 - PCA transform")
scatter(X_PCA,title="X_PCA - PCA reconstruct")

scatter(Y2,title="Y2 - PPCA transform")
pppca = scatter(X_PPCA,title="X_PPCA - PPCA reconstruct");
plot(pdata,pppca)

#####################################################################################################

using ManifoldLearning
Y_Isomap = transform(Isomap,data_array; k=12,d=2)

Y_DiffMap = transform(DiffMap, data_array; d=2, t=1, ɛ=1.0)

Y_LEM = transform(LEM,data_array; k=12, d=2, t=1.0)

Y_LLE = transform(LLE, data_array; k = 12, d = 2)

Y_HLLE = transform(HLLE, data_array; k = 12, d = 2)

Y_LTSA = transform(LTSA, data_array; k = 12, d = 2)

print(Y_LEM)
plot(Y_LTSA)
scatter(Y_LEM, title="LEM - ManifoldLearning")
