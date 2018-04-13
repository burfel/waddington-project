using MultivariateStats

# suppose Xtr and Xte are training and testing data matrix,
Xtr = rand(200, 300)
Xte = rand(200, 300)
# with each observation in a column

# train a PCA model
M = fit(PCA, Xtr; maxoutdim=100)

# apply PCA model to testing set
Yte = transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = reconstruct(M, Yte)

using RDatasets
using Plots
using DataArrays


plotly() # using plotly for 3D-interacive graphing

# load iris dataset

iris = dataset("datasets", "iris")

# split half to training set
# Xtr = convert(Array,DataArrays.DataArray(iris[1:2:end,1:4]))'
Xtr = convert(Array,DataFrames.DataFrame(iris[1:2:end,1:4]))'
# Xtr_labels = convert(Array,DataArray(iris[1:2:end,5]))
Xtr_labels = convert(Array,DataFrames.DataFrame(iris[1:2:end,5]))

# split other half to testing set
# Xte = convert(Array,DataArray(iris[2:2:end,1:4]))'
Xte = convert(Array,DataFrames.DataFrame(iris[2:2:end,1:4]))'
Xte_labels = convert(Array,DataFrames.DataFrame(iris[2:2:end,5]))

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a PCA model, allowing up to 3 dimensions
M = fit(PCA, Xtr; maxoutdim=3)

# apply PCA model to testing set
Yte = transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = reconstruct(M, Yte)

##############---runs till here
# group results by testing set labels for color coding
rel = Xte_labels.=="setosa"
setosa = Yte[names(Yte)[rel]
"""
setosa = Yte[:,Xte_labels.=="setosa"]
versicolor = Yte[:,Xte_labels.=="versicolor"]
virginica = Yte[:,Xte_labels.=="virginica"]
"""

# visualize first 3 principal components in 3D interacive plot
p = scatter(setosa[1,:],setosa[2,:],setosa[3,:],marker=:circle,linewidth=0)
scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")
