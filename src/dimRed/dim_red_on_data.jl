# -*- coding: utf-8 -*-
#=
The following program tries to implement the dimensionality reduction functions
available in the packages MultivariateStats.jl and ManifoldLearning.jl

Note: You have to run Read.jl first (dataframe data needed).

Args:
    -

Returns:
    -
=#

using DataStructures
using MultivariateStats
using ManifoldLearning
using Plots
using PyPlot


Base.compilecache("MultivariateStats")
Pkg.checkout("MultivariateStats")


data_array = convert(Array, data)
#data_array = collect(skipmissing(data_array))
data_array = reshape(data_array, (547,96))


# ---- DIMENSIONALITY REDUCTION using MultivariateStats.jl-----------------------

#----- 1. PRINCIPAL COMPONENT ANALYSIS (PCA)-------------------------------------

# Observations have to equal columns, ie we need to transpose data_array
data_arrayT = transpose(data_array)

# Perform PCA over the data given in a matrix data_array. Each column of data_array is an observation.
# Returns an instance of PCA.
M1_max = fit(PCA, data_arrayT; maxoutdim=96)
M1 = fit(PCA, data_arrayT; maxoutdim=2)

# transforms observations data_array into PCs
Y1_max = transform(M1_max, data_arrayT)
Y1 = transform(M1, data_arrayT)

X_PCA_max = reconstruct(M1_max, Y1_max)
X_PCA = reconstruct(M1, Y1)

# Get the projection matrix (of size (d, p)).
# Each column of the projection matrix corresponds to a principal component.
# The principal components are arranged in descending order of the corresponding variances.
M1_proj_max = projection(M1_max)
M1_proj = projection(M1)
#indim(M1_max)
#outdim(M1_max)
#indim(M1)
#outdim(M1)

# Get the input dimension d, i.e the dimension of the observation space.
indim(M1_max)
indim(M1)

# Get the output dimension p, i.e the dimension of the principal subspace.
outdim(M1_max) #75
outdim(M1)

# The variances of principal components.
pca_var_max = principalvars(M1_max)
Plots.bar(pca_var_max, xlabel = "principal components")
PyPlot.title("Eigenvalues / Explained variance of principal components")
PyPlot.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs_max")

pca_var = principalvars(M1)
Plots.bar(pca_var, xlabel = "principal components")
PyPlot.title("Eigenvalues / Explained variance of principal components")
PyPlot.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs")

# The total variance of principal components, which is equal to sum(principalvars(M)).
print("total variance of PCs: ", tprincipalvar(M1_max)) # 1643.3876618996323
print("total variance of PCs: ", tprincipalvar(M1)) # 501.42985340006794

# The total residual variance.
print("total residual variance: ", tresidualvar(M1_max)) # 15.604526092502056
print("total residual variance: ", tprincipalvar(M1)) # 501.42985340006794

# The total observation variance, which is equal to tprincipalvar(M) + tresidualvar(M).
print("total observation variance: ", tvar(M1_max)) # 1658.9921879921344
print("total observation variance: ", tvar(M1)) # 1658.9921879921344

# The ratio of variance preserved in the principal subspace, which is equal to tprincipalvar(M) / tvar(M).
print("ratio of variance preserved in the principal subspace: ", principalratio(M1_max)) # 0.9905939725301611
print("ratio of variance preserved in the principal subspace: ",principalratio(M1)) # 0.3022496772615576


# PLOTS
Plots.scatter(data_array, legend=false)
#PyPlot.title("Data set")
#PyPlot.savefig("../Single_cell_data/plots/data_set")

# PCA transform
Plots.scatter(Y1_max,title="PCA transform", legend=false)
PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_transform_max")

Plots.scatter(Y1,title="PCA transform", legend=false) # compare it to own function
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_transform")

Plots.scatter(X_PCA_max,title="PCA reconstruct", legend=false)# compare it to own function
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct_max")

Plots.scatter(X_PCA,title="PCA reconstruct", legend=false)
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct")

#=
# projected onto 75 PCs
Plots.scatter(M1_proj_max,title="PCs")
#PyPlot.title# projected onto 75 PCs
Plots.scatter(M1_proj_max,title="PCs")
#PyPlot.title("PCA reconstruct")
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct")

# projected onto 2 PCs
Plots.scatter(M1_proj,title="PCs")("PCA reconstruct")
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct")

# projected onto 2 PCs
Plots.scatter(M1_proj,title="PCs")
=#

####-------------

#= not running:
using LowRankModels
import ScikitLearnBase

ScikitLearnBase.fit_transform!(LowRankModels.PCA(k=3, max_iter=500), data_array)
=#

#----- 2. Classical Multidimensional Scaling (MDS)-------------------------------
# This function derives a p-dimensional embedding based on a given distance_matrix.
# It returns a coordinate matrix of size (p, n), where each column is the coordinates for an observation.

# gram matrix  = X * X_transpose
gram_matrix = data_array * transpose(data_array)
# compute the distance matrix
distance_matrix = gram2dmat(gram_matrix)
MDS = classical_mds(distance_matrix, 96, dowarn=true)
MDSt = transpose(MDS)

Plots.scatter(MDSt,title="Multidimensional scaling (MDS)", legend=false)
PyPlot.savefig("../Single_cell_data/plots/MDS/mds")


#----- 3. Linear Discriminant Analysis (LDA)-------------------------------------
# Multiclass LDA (Linear Discriminant Analysis)
#=
Multi-class LDA is a generalization of standard two-class LDA that can handle arbitrary number of classes.
Linear Discriminant Analysis are statistical analysis methods to find a linear combination of features for
separating observations in n classes.
=#
#fit(MulticlassLDA, number_of_classes, data_array, vector_of_class_labels; ...)
#Todo: define number_of_classes, vector_of_class_labels


#----- 3. Independent Component Analysis (ICA)-------------------------------------
# Independent Component Analysis (ICA), FastICA
#=
Independent Component Analysis (ICA) is a computational technique for separating a multivariate signal
into additive subcomponents, with the assumption that the subcomponents are non-Gaussian and independent
from each other.
----DOES NOT RUN-----FAILS TO CONVERGE after 100 iterations-----------------------
=#
M12 = fit(ICA, data_array, 2)
Y12 = transform(M12, data_array)
X_ICA = reconstruct(M12, Y12)

#scatter(Y12, title="ICA transform")
#scatter(X_ICA, title="ICA reconstruct")


#----- 4. Probabilistic PCA (PPCA)-----------------------------------------------
M2 = fit(PPCA, data_array; maxoutdim=96)
Y2 = transform(M2, data_array)
#PrincComp2 = projection(data_array)
X_PPCA = reconstruct(M2, Y2)

# PLOTS
Plots.scatter(Y2,title="PPCA transform", legend=false)
#PyPlot.title("PPCA transform")
#PyPlot.savefig("../Single_cell_data/plots/PPCA/PPCA_transform")

Plots.scatter(X_PPCA,title="PPCA reconstruct", legend=false)
#PyPlot.title("PPCA reconstruct")
#PyPlot.savefig("../Single_cell_data/plots/PPCA/PPCA_reconstruct")


#----- 5. Factor Analysis -------------------------------------------------------
#=Factor Analysis (FA) is a linear-Gaussian latent variable model that is closely
related to probabilistic PCA. In contrast to the probabilistic PCA model, the covariance
of conditional distribution of the observed variable given the latent variable is
diagonal rather than isotropic.
----DOES NOT RUN----- DOES NOT TERMINATE (within a considerable amount of time)
=#
# train a FactorAnalysis model
M22 = fit(FactorAnalysis, data_array; maxoutdim=2)

# apply FactorAnalysis model to testing set
Y22 = transform(M22, data_array)

# reconstruct testing observations (approximately)
X_FactorAna = reconstruct(M22, Y22)


#----- 6. Kernel PCA ------------------------------------------------------------

# train a kernel PCA model
M3_max = fit(KernelPCA, data_array; maxoutdim=96, inverse=true)
M3 = fit(KernelPCA, data_array; maxoutdim=2, inverse=true)

# apply kernel PCA model to testing set
Y3_max = transform(M3_max, data_array)
Y3 = transform(M3, data_array)

# reconstruct testing observations (approximately)
X_kernel_max = reconstruct(M3_max, Y3_max)
X_kernel = reconstruct(M3, Y3)

# PLOTS
Plots.scatter(Y3_max, title="kernelPCA transform (max)", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCAm_transform")

Plots.scatter(X_kernel_max, title="kernelPCA reconstruct (max)", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCAm_reconstruct")

Plots.scatter(Y3, title="kernelPCA transform", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCA_transform")

Plots.scatter(X_kernel, title="kernelPCA reconstruct", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCA_reconstruct")


###--------------------------------------------------------------------------------

# ---- DIMENSIONALITY REDUCTION using ManifoldLearning.jl-------------------------

# Isomap
Y_Isomap = transform(Isomap,data_array; k=12,d=2)

# Diffusion maps
Y_DiffMap = transform(DiffMap, data_array; d=2, t=1, ɛ=1.0)

# Laplacian Eigenmaps -- working
Y_LEM = transform(LEM,data_array; k=12, d=2, t=1.0)
Plots.scatter(Y_LEM, title="ManifoldLearning: Laplacian Eigenmaps")
#PyPlot.savefig("../Single_cell_data/plots/manifold-learning/lem")

# Local linear embedding
Y_LLE = transform(LLE, data_array; k = 12, d = 2)

# Hessian Eigenmaps transformation
Y_HLLE = transform(HLLE, data_array; k = 12, d = 2)

# Local tangent space alignment -- working
Y_LTSA = transform(LTSA, data_array; k = 12, d = 2)
Plots.scatter(Y_LTSA, title="ManifoldLearning: Local tangent space alignment")
#PyPlot.savefig("../Single_cell_data/plots/manifold-learning/ltsa")


#--------other useful methods
# Canonical Correlation Analysis (CCA)
#=
Canonical Correlation Analysis (CCA) is a statistical analysis technique to
identify correlations between two sets of variables. Given two vector variables X
and Y, it finds two projections, one for each, to transform them to a common space
with maximum correlations.
=#
