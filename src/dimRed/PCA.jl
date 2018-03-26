using MultivariateStats, RDatasets, Plots
#plotly() # using plotly for 3D-interacive graphing

#=
# load iris dataset
iris = dataset("datasets", "iris")

# split half to training set
Xtr = convert(Array,DataArray(iris[1:2:end,1:4]))'
Xtr_labels = convert(Array,DataArray(iris[1:2:end,5]))

# split other half to testing set
Xte = convert(Array,DataArray(iris[2:2:end,1:4]))'
Xte_labels = convert(Array,DataArray(iris[2:2:end,5]))

# suppose Xtr and Xte are training and testing data matrix,
# with each observation in a column

# train a PCA model, allowing up to 3 dimensions
M = fit(PCA, Xtr; maxoutdim=3)

# apply PCA model to testing set
Yte = transform(M, Xte)

# reconstruct testing observations (approximately)
Xr = reconstruct(M, Yte)

# group results by testing set labels for color coding
setosa = Yte[:,Xte_labels.=="setosa"]
versicolor = Yte[:,Xte_labels.=="versicolor"]
virginica = Yte[:,Xte_labels.=="virginica"]

# visualize first 3 principal components in 3D interacive plot
p = Plots.scatter(setosa[1,:],setosa[2,:],setosa[3,:],marker=:circle,linewidth=0)
scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")
=#

#Xtr_labels = convert(Array,DataArray(iris[1:2:end,5]))
#M1_max = fit(PCA, data_arrayT; maxoutdim=96)
M1 = fit(PCA, data_arrayT; maxoutdim=3)

# transforms observations data_array into PCs
#Y1_max = transform(M1_max, data_arrayT) # 75x547
Y1 = transform(M1, data_arrayT) # 2x547

#X_PCA_max = reconstruct(M1_max, Y1_max) # 96x547
X_PCA = reconstruct(M1, Y1) # 96x547

# Get the projection matrix (of size (d, p)).
# Each column of the projection matrix corresponds to a principal component.
# The principal components are arranged in descending order of the corresponding variances.
#M1_proj_max = projection(M1_max) # 96x75
M1_proj = projection(M1) # 96x2


##-------------PLOT TRANSFORMATION (projection onto PCs)-------------------------------
# PROJECTION ONTO 3 PCs
p = Plots.scatter(Y1[1,:],Y1[2,:],Y1[3,:],marker=:circle,linewidth=0, title="PCA projection onto the first 3 PCs")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(p,xlabel="PC1 (21.85%)",ylabel="PC2 (8.38%)",zlabel="PC3 (6.56%)", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.pcasavefig("PCA_3D")

# PROJECTION ONTO 2 PCs
q = Plots.scatter(Y1[1,:],Y1[2,:],marker=:circle,linewidth=0, title="PCA projection onto the first 2 PCs")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(q,xlabel="PC1 (21.85%)",ylabel="PC2 (8.38%)", legend=false)
#PyPlot.title("PCA projection onto the first 2 PCs")
Plots.savefig("PCA_2D")

##--------------STATISTICS-----------------------------------------------------------------
principalvars(M1)

principalratio(M1)
tprincipalvar(M1)
tresidualvar(M1)

# principal ratio preserved in the subspace: 36.8
# 1st PC: 21.85%
# 2nd PC: 8.38%
# 3rd PC: 6.56%

##-------------PLOT RECONSTRUCTION--------------------------------------------------
r = Plots.scatter(X_PCA[31,:],X_PCA[22,:],X_PCA[11,:],marker=:circle,linewidth=0, title="PCA reconstruction from the first 3 PCs")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(r,xlabel="Gdf3",ylabel="Fgf4",zlabel="Cldn6", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("PCA_reconstruction_3D_impGenes")

s = Plots.scatter(X_PCA[31,:],X_PCA[22,:], marker=:circle,linewidth=0, title="PCA reconstruction from the first 2 PCs")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(s,xlabel="Gdf3",ylabel="Fgf4", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("PCA_reconstruction_2D_impGenes")

##-------------old data set for comparison-------------------------------------
uu = Plots.scatter(data_array[31,:], data_array[22,:],data_array[11,:],marker=:circle,linewidth=0, title="Visualisation of the original data set (3 dimensions with highest variance)")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(uu,xlabel="Gdf3",ylabel="Fgf4",zlabel="Cldn6", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("PCA_reconstruction_3D_impGenes")

u = Plots.scatter(data_array[31,:],data_array[22,:], marker=:circle,linewidth=0, title="Visualisation of the original data set (2 dimensions with highest variance)")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(u,xlabel="Gdf3",ylabel="Fgf4", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
#Plots.savefig("PCA_reconstruction_2D_impGenes")


vv = Plots.scatter(data_array[31,:], data_array[22,:],data_array[11,:],marker=:circle,linewidth=0, title="Visualisation of the original data set (first 3 dimensions)")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(vv,xlabel="x_1",ylabel="x_2",zlabel="x_3", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
#Plots.savefig("PCA_reconstruction_3D_impGenes")

v = Plots.scatter(data_array[31,:],data_array[22,:], marker=:circle,linewidth=0, title="Visualisation of the original data set (first 2 dimensions)")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(v,xlabel="x_1",ylabel="x_2", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
#Plots.savefig("PCA_reconstruction_2D_impGenes")


##-------------PLOT PROJECTION (rotation matrix)---------------------------------------
t = Plots.bar(M1_proj[1,:],M1_proj[2,:],M1_proj[3,:],marker=:circle,linewidth=0, title="PCA rotation matrix for the first 3 PCs")
#scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
#scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
Plots.plot!(t,xlabel="x_1",ylabel="x_2",zlabel="x_3", legend=false)
#PyPlot.title("PCA projection onto the first 3 PCs")
Plots.savefig("PCA_rotation_3D")
