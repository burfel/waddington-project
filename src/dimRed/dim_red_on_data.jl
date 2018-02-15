#=
The following program tries to implement the dimensionality reduction functions
available in the packages MultivariateStats.jl and ManifoldLearning.jl
Note: You have to run topCorr.jl first (dataframe data needed).
=#
using DataStructures
using MultivariateStats
using ManifoldLearning
using Plots

data_array = convert(Array, data)
#data_array = collect(skipmissing(data_array))
data_array = reshape(data_array, (547,96))


# ----MultivariateStats.jl DIMENSIONALITY REDUCTION

# Principal Components Analysis (PCA)
M1 = fit(PCA, data_array; maxoutdim=2)
Y1 = transform(M1, data_array)
X_PCA = reconstruct(M1, Y1)
projection(X_PCA)

#= not trunning:
using LowRankModels
import ScikitLearnBase

ScikitLearnBase.fit_transform!(LowRankModels.PCA(k=3, max_iter=500), data_array)
=#

# Classical Multidimensional Scaling (MDS)
#distance_matrix = gram2dmat(gram_matrix)
#classical_mds(distance_matrix, 2[, dowarn=true])
#Todo: compute gram_matrix

# Linear Discriminant Analysis (LDA)

# Multiclass LDA (Linear Discriminant Analysis)
#=
Multi-class LDA is a generalization of standard two-class LDA that can handle arbitrary number of classes.
Linear Discriminant Analysis are statistical analysis methods to find a linear combination of features for
separating observations in n classes.
=#
#fit(MulticlassLDA, number_of_classes, data_array, vector_of_class_labels; ...)
#Todo: define number_of_classes, vector_of_class_labels

# Independent Component Analysis (ICA), FastICA
#=
Independent Component Analysis (ICA) is a computational technique for separating a multivariate signal
into additive subcomponents, with the assumption that the subcomponents are non-Gaussian and independent
from each other.
=#
M12 = fit(ICA, data_array, 2)
Y12 = transform(M12, data_array)
X_ICA = reconstruct(M12, Y12)

# Probabilistic PCA
M2 = fit(PPCA, data_array; maxoutdim=2)
Y2 = transform(M2, data_array)
#PrincComp2 = projection(data_array)
X_PPCA = reconstruct(M2, Y2)

# Factor Analysis
#=Factor Analysis (FA) is a linear-Gaussian latent variable model that is closely
related to probabilistic PCA. In contrast to the probabilistic PCA model, the covariance
of conditional distribution of the observed variable given the latent variable is
diagonal rather than isotropic.
=#
# train a FactorAnalysis model
M22 = fit(FactorAnalysis, Xtr; maxoutdim=100)
# apply FactorAnalysis model to testing set
Y22 = transform(M22, data_array)
# reconstruct testing observations (approximately)
X_FactorAna = reconstruct(M22, Y22)

# Kernel PCA
# train a kernel PCA model
M3 = fit(KernelPCA, data_array; maxoutdim=100, inverse=true)
# apply kernel PCA model to testing set
Y3 = transform(M3, data_array)
# reconstruct testing observations (approximately)
X_kernel = reconstruct(M3, Y3)


# PLOTTING linear methods
scatter(data_array)

scatter(Y1,title="Y1 - PCA transform")
scatter(X_PCA,title="X_PCA - PCA reconstruct")

scatter(Y2,title="Y2 - PPCA transform")
scatter(X_PPCA,title="X_PPCA - PPCA reconstruct")

scatter(Y12,title="Y12 - ICA transform")
#scatter(X_ICA,title="X_ICA - ICA reconstruct")


# ----ManifoldLearnin.jl DIMENSIONALITY REDUCTION

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


# PLOTTING nonlinear methods
scatter(Y_LEM, title="LEM - ManifoldLearning")

scatter(Y_LTSA, title="LTSA - ManifoldLearning")


#--------other useful methods
# Canonical Correlation Analysis (CCA)
#=
Canonical Correlation Analysis (CCA) is a statistical analysis technique to
identify correlations between two sets of variables. Given two vector variables X
and Y, it finds two projections, one for each, to transform them to a common space
with maximum correlations.
=#
