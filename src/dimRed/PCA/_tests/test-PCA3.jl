using DataFrames

# Some example data
data2 = convert(DataFrame, [10  10  10  8   8.04   9.14  7.46   6.58;
 8   8   8   8   6.95   8.14  6.77   5.76;
 13  13  13  8   7.58   8.74  12.74  7.71;
 9   9   9   8   8.81   8.77  7.11   8.84;
 11  11  11  8   8.33   9.26  7.81   8.47;
 14  14  14  8   9.96   8.1   8.84   7.04;
 6   6   6   8   7.24   6.13  6.08   5.25;
 4   4   4   19  4.26   3.1   5.39   12.5;
 12  12  12  8   10.84  9.13  8.15   5.56;
 7   7   7   8   4.82   7.26  6.42   7.91;
 5   5   5   8   5.68   4.74  5.73   6.89]);

# This returns the column headings as an array of symbols
n = names(data2)

# Here is some made up Boolean array, similar to what you had
present = [true, true, false, false, true, true, true, false]

# Then just filter the dataframe to only return the correct columns
data2[names(data2)[present]]


#-----------

"""
An Array{Bool} stores each true/false value as a Bool, which is represented internally as a UInt8. So if your array has N elements, it will take N bytes to store it.

A BitArray stores each true/false value as a single bit, with (conceptually) 8 of them packed into a single UInt8. Consequently, it takes only N/8 bytes to store the array. A BitArray also has methods defined that handle all the required bit-twiddling operations for you.
"""

# linear DR algorithms

##PCA

using MultivariateStats

# suppose Xtr and Xte are training and testing data matrix,# with each observation in a column
# train a PCA modelM = fit(PCA, Xtr; maxoutdim=100)
# apply PCA model to testing setYte = transform(M, Xte)
# reconstruct testing observations (approximately)Xr = reconstruct(M, Yte)

using MultivariateStats, RDatasets, Plots
plotly() # using plotly for 3D-interacive graphing

# load iris dataset
iris = dataset("datasets", "iris")
data = (DataFrame, iris)

# split half to training set
Xtr = convert(Array,DataFrames.DataFrame(data[1:2:end,1:4]))'
Xtr_labels = convert(Array,DataFrames.DataFrame(data[1:2:end,5]))

# split other half to testing set
Xte = convert(Array,DataFrames.DataFrame(data[2:2:end,1:4]))'
Xte_labels = convert(Array,DataFrames.DataFrame(data[2:2:end,5]))

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a PCA model, allowing up to 3 dimensions
M = MultivariateStats.fit(PCA, Xtr; maxoutdim=3)

# apply PCA model to testing set
Yte = MultivariateStats.transform(M, Xte)
# reconstruct testing observations (approximately)
Xr = MultivariateStats.reconstruct(M, Yte)

# group results by testing set labels for color coding
data = convert(DataFrame, Yte)
rel_setosa = transpose(Xte_labels.=="setosa")
setosa = data[names(data)[rel_setosa[:,1]]]

rel_versicolor = transpose(Xte_labels.=="versicolor")
versicolor = data[names(data)[rel_versicolor[:,1]]]

rel_virginica = transpose(Xte_labels.=="virginica")
virginica = data[names(data)[rel_virginica[:,1]]]

#versicolor = Yte[:,Xte_labels.=="versicolor"[:,1]]
#virginica = Yte[:,Xte_labels.=="virginica"]


#####--- runs till here
"""
# visualize first 3 principal components in 3D interacive plot
p = scatter(setosa[1,:],setosa[2,:],setosa[3,:],marker=:circle,linewidth=0)
scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")

@df setosa scatter([1,:], [2])
"""
"""
###PPCA
#using MultivariateStats

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a PCA model
M = MultivariateStats.fit(PPCA, Xtr; maxoutdim=100)

# apply PCA model to testing set
Yte = MultivariateStats.transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = MultivariateStats.reconstruct(M, Yte)

#--------

###kPCA
#using MultivariateStats

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a kernel PCA model
#M = MultivariateStats.fit(KernelPCA, Xtr; maxoutdim=100, inverse=true)

# apply kernel PCA model to testing set
Yte = MultivariateStats.transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = MultivariateStats.reconstruct(M, Yte)

#---------

# FA

using MultivariateStats

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a FactorAnalysis model
M = MultivariateStats.fit(FactorAnalysis, Xtr; maxoutdim=100)

# apply FactorAnalysis model to testing set
Yte = MultivariateStats.transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = MultivariateStats.reconstruct(M, Yte)


#----MANIFOLD LEARNING

# Isomap
using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply Isomap transformation to the dataset
Y = ManifoldLearning.transform(Isomap, X; k = 12, d = 2)

# Diffsion map
#using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply DiffMap transformation to the dataset
Y = ManifoldLearning.transform(DiffMap, X; d=2, t=1, ɛ=1.0)


# Laplacian Eigenmaps
#using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply Laplacian Eigenmaps transformation to the dataset
Y = ManifoldLearning.transform(LEM, X; k = 12, d = 2, t = 1.0)


# HLLE
#using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply LLE transformation to the dataset
Y = ManifoldLearning.transform(LLE, X; k = 12, d = 2)


# HHE
#using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply HLLE transformation to the dataset
Y = ManifoldLearning.transform(HLLE, X; k = 12, d = 2)


# LTSA
#using ManifoldLearning

# suppose X is a data matrix, with each observation in a column
# apply LTSA transformation to the dataset
Y = ManifoldLearning.transform(LTSA, X; k = 12, d = 2)
"""
"""
potential things to do:
- ICA with other algorithm
- topslam in python --> call in Julia
- embed skikit functions in julia
- (BGLPVM too much)
"""
