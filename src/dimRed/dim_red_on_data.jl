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
#using PyPlot


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
M1_3 = fit(PCA, data_arrayT; maxoutdim=3)

# transforms observations data_array into PCs
Y1_max = transform(M1_max, data_arrayT) # 75x547
Y1 = transform(M1, data_arrayT) # 2x547
Y1_3 = transform(M1_3, data_arrayT)

X_PCA_max = reconstruct(M1_max, Y1_max) # 96x547
X_PCA = reconstruct(M1, Y1) # 96x547
X_PCA = reconstruct(M1_3, Y1_3)

# Get the projection matrix (of size (d, p)).
# Each column of the projection matrix corresponds to a principal component.
# The principal components are arranged in descending order of the corresponding variances.
M1_proj_max = projection(M1_max) # 96x75
M1_proj = projection(M1) # 96x2
M1_proj = projection(M1_3) # 96x3
#indim(M1_max)
#outdim(M1_max)
#indim(M1)
#outdim(M1)

##---LINEAR COMBINATIONS OF OLD BASIS (GENES)
pc1 = M1_proj[:,1]
Plots.bar(pc1, ylabel="weight of particular gene", xlabel="genes (1-96)", legend=false, title="Composition of first PC")
Plots.savefig("../Single_cell_data/plots/PCA/compPC1")

pc2 = M1_proj[:,2]
Plots.bar(pc2, ylabel="weight of particular gene", xlabel="genes (1-96)", legend=false, title="Composition of second PC")
Plots.savefig("../Single_cell_data/plots/PCA/compPC2")

pc3 = M1_proj[:,3]
Plots.bar(pc3, ylabel="weight of particular gene", xlabel="genes (1-96)", legend=false, title="Composition of third PC")
Plots.savefig("../Single_cell_data/plots/PCA/compPC3")
##----
# WHICH ARE THE BIGGEST WEIGHTS?
names = DataFrames.names(data)
#max_pc1=findmax(pc1)
#names[71] #Rai1
min_pc1=findmin(pc1) #
names[22] #Fgf4

pc_total = sum(transpose(abs(M1_proj_max)),1)
Plots.bar(transpose(pc_total), ylabel="weight of particular gene", xlabel="genes (1-96)", legend=false, title="Contribution of genes summed over all PCs")
Plots.savefig("../Single_cell_data/plots/PCA/compPC_total")

max_pc2=findmax(pc2)
names[11] #Cldn6
#min_pc2=findmin(pc2)

#max_pc3=findmax(pc3)
#names[45] #Klf4
min_pc3=findmin(pc3)
names[8] #Cdh2

print(data[:Smarca4])

findmax(pc_total)
findmin(pc_total)
names[89] #Tubb3
names[49] #MBP

# Get the input dimension d, i.e the dimension of the observation space.
indim(M1_max) #96
indim(M1) #96
indim(M1_3) #96

# Get the output dimension p, i.e the dimension of the principal subspace.
outdim(M1_max) #75
outdim(M1) #2
outdim(M1_3) #3

# The variances of principal components.
pca_var_max = principalvars(M1_max)
Plots.bar(pca_var_max, xlabel = "Principal components", ylabel = "Eigenvalues", legend=false, title="Explained variances of all the PCs")
Plots.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs_max")

# NORMALISED variance:

maxVar = tvar(M1_max)
normalise(x) = x/maxVar
normalisedEVs = normalise.(pca_var_max)
print(normalisedEVs)

cum=cumsum(normalisedEVs)
Plots.bar(normalisedEVs, xlabel = "Principal components", ylabel = "Explained variance", legend=false, title="Explained variances of all the PCs")
Plots.plot!(cum, label="cumulative variance")
Plots.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs_max_normed_cum")
#savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs_max_normed")

pca_var = principalvars(M1)
Plots.bar(pca_var, xlabel = "Principal components", ylabel = "Eigenvalues", legend=false, title="Explained variances of 2 PCs")
Plots.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs")

pca_var_3 = principalvars(M1_3)
Plots.bar(pca_var_3, xlabel = "Principal components", ylabel = "Eigenvalues", legend=false, title="Explained variances of 3 PCs")
Plots.savefig("../Single_cell_data/plots/PCA/eigenvalues_of_resp_PCs_3")

# The total variance of principal components, which is equal to sum(principalvars(M)).
print("total variance of PCs: ", tprincipalvar(M1_max)) # 1643.3876618996323
print("total variance of PCs: ", tprincipalvar(M1)) # 501.42985340006794
print("total variance of PCs: ", tprincipalvar(M1_3)) # 610.3081101099062

# The total residual variance.
print("total residual variance: ", tresidualvar(M1_max)) # 15.604526092502056
print("total residual variance: ", tprincipalvar(M1)) # 501.42985340006794
print("total residual variance: ", tprincipalvar(M1_3)) # 610.3081101099062

# The total observation variance, which is equal to tprincipalvar(M) + tresidualvar(M).
print("total observation variance: ", tvar(M1_max)) # 1658.9921879921344
print("total observation variance: ", tvar(M1)) # 1658.9921879921344
print("total observation variance: ", tvar(M1_3)) # 1658.9921879921344

# The ratio of variance preserved in the principal subspace, which is equal to tprincipalvar(M) / tvar(M).
print("ratio of variance preserved in the principal subspace: ", principalratio(M1_max)) # 0.9905939725301611
print("ratio of variance preserved in the principal subspace: ",principalratio(M1)) # 0.3022496772615576
print("ratio of variance preserved in the principal subspace: ",principalratio(M1_3)) # 0.36787883302124375


# PLOTS (not very meaningful)
Plots.scatter(data_array[:,1:5], legend=false, title="Sample data set", xlabel="Samples (1-547)", ylabel="Expression level of genes 1-5", legend=true)
Plots.plot(data_array, legend=false)
Plots.savefig("../Single_cell_data/plots/data_set_5genes")
#PyPlot.savefig("../Single_cell_data/plots/data_set")

Plots.scatter(transpose(data_array), legend=false, title="Sample data set", xlabel="Genes (1-96)", ylabel="Gene expression level", legend=false)
Plots.plot(data_array, legend=false)
Plots.savefig("../Single_cell_data/plots/data_set_T")
#Plots.savefig("../Single_cell_data/plots/data_set")
#PyPlot.savefig("../Single_cell_data/plots/data_set")

#=
# PCA transform
Plots.scatter(Y1_max,title="PCA transform", legend=false)
# Plots.plot(Y1_max,title="PCA transform", legend=false)
PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_transform_max")

Plots.scatter(Y1,title="PCA transform", legend=false) # compare it to own function
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_transform")
Plots.plot(Y1,title="PCA transform", legend=false)

Plots.scatter(M1_proj_max,title="PCA projection", legend=false) # compare it to own function
PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_projection_max")
#Plots.plot(M1_proj_max,title="PCA projection", legend=false)

Plots.scatter(M1_proj,title="PCA projection", legend=false) # compare it to own function
PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_projection")

Plots.scatter(X_PCA_max,title="PCA reconstruct", legend=false)# compare it to own function
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct_max")

Plots.scatter(X_PCA,title="PCA reconstruct", legend=false)
#PyPlot.savefig("../Single_cell_data/plots/PCA/PCA_reconstruct")

=#

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

#= PLOT (not very meaningful)
Plots.scatter(MDSt,title="Multidimensional scaling (MDS)", legend=false)
savefig("../Single_cell_data/plots/MDS/mds")
=#

# PROJECTION ONTO 3 new directions
z = Plots.scatter(MDS[1,:],MDS[2,:],MDS[3,:],marker=:circle,linewidth=0, title="Classical Multidimensional scaling")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(z,xlabel="x_new",ylabel="y_new",zlabel="z_new", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
savefig("../Single_cell_data/plots/MDS/mds_3D")



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

#=
# PLOTS (not very meaningful)
Plots.scatter(Y2,title="PPCA transform", legend=false)
#PyPlot.title("PPCA transform")
#PyPlot.savefig("../Single_cell_data/plots/PPCA/PPCA_transform")

Plots.scatter(X_PPCA,title="PPCA reconstruct", legend=false)
#PyPlot.title("PPCA reconstruct")
#PyPlot.savefig("../Single_cell_data/plots/PPCA/PPCA_reconstruct")
=#

# PLOT PROJECTION
x = Plots.scatter(Y2[1,:],Y2[2,:], marker=:circle,linewidth=0, title="PPCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(x,xlabel="x_new",ylabel="y_new", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/PPCA/ppca_2D")

# PLOT RECONSTRUCTION 3D
a = Plots.scatter(X_PPCA[1,:],X_PPCA[2,:], X_PPCA[3,:], marker=:circle,linewidth=0, title="PPCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(a,xlabel="x_1",ylabel="x_2", zlabel="x_3", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/PPCA/ppca_3D_recon")

# PLOT RECONSTRUCTION 2D
b = Plots.scatter(X_PPCA[1,:],X_PPCA[2,:], marker=:circle,linewidth=0, title="PPCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(b,xlabel="x_1",ylabel="x_2", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/PPCA/ppca_2D_recon")


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
#M3_max = fit(KernelPCA, data_array; maxoutdim=96, inverse=true)
M3 = fit(KernelPCA, data_array; maxoutdim=2, inverse=true)

# apply kernel PCA model to testing set
#Y3_max = transform(M3_max, data_array)
Y3 = transform(M3, data_array)

# reconstruct testing observations (approximately)
#X_kernel_max = reconstruct(M3_max, Y3_max)
X_kernel = reconstruct(M3, Y3)

#=
# PLOTS (not very meaningful)
Plots.scatter(Y3_max, title="kernelPCA transform (max)", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCAm_transform")

Plots.scatter(X_kernel_max, title="kernelPCA reconstruct (max)", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCAm_reconstruct")

Plots.scatter(Y3, title="kernelPCA transform", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCA_transform")

Plots.scatter(X_kernel, title="kernelPCA reconstruct", legend=false)
PyPlot.savefig("../Single_cell_data/plots/kernelPCA/kernelPCA_reconstruct")
=#

c = Plots.scatter(Y3[1,:],Y3[2,:], marker=:circle,linewidth=0, title="Kernel PCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(c,xlabel="PC1",ylabel="PC2", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/kernelPCA/kpca_2D")


d = Plots.scatter(X_kernel[1,:],X_kernel[2,:], marker=:circle,linewidth=0, title="Reconstruction from Kernel PCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(d,xlabel="x_1",ylabel="x_2", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/kernelPCA/kpca_recon_2D")


e = Plots.scatter(X_kernel[1,:],X_kernel[2,:], X_kernel[3,:],marker=:circle,linewidth=0, title="Reconstruction from Kernel PCA")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(e,xlabel="x_1",ylabel="x_2", zlabel="x_3", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("../Single_cell_data/plots/kernelPCA/kpca_recon_3D")



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


#= ----cannot plot object LEM (no user recipe defined)----------
f = Plots.scatter(Y_LEM[1,:],Y_LEM[2,:], Y_LEM[3,:],marker=:circle,linewidth=0, title="Laplacian Eigenmaps")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(f,xlabel="x_1",ylabel="x_2", zlabel="x_3", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("LEM")
=#


#--------other useful methods
# Canonical Correlation Analysis (CCA)
#=
Canonical Correlation Analysis (CCA) is a statistical analysis technique to
identify correlations between two sets of variables. Given two vector variables X
and Y, it finds two projections, one for each, to transform them to a common space
with maximum correlations.
=#
