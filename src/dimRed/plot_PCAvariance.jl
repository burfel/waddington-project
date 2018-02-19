#=
The following program plots the explained variance of the data points in data_array
Note: You have to run topCorr.jl first (dataframe data needed).
=#

using ScikitLearn
using Plots
using PyPlot
using ScikitLearn.GridSearch: GridSearchCV
using ScikitLearn.Pipelines: Pipeline, named_steps

@sk_import decomposition: PCA
@sk_import datasets: load_digits
@sk_import linear_model: LogisticRegression

# Plot the PCA spectrum
fit!(pca, data_array)

figure(1, figsize=(4, 3))
clf()
axes([.2, .2, .7, .7])
plot(pca[:explained_variance_], linewidth=2)
axis("tight")
xlabel("n components")
ylabel("explained variance")


# SVD without ScikitLearn:
F = svdfact(data_array)
#F[:U]
#F[:S]
#F[:V]
singular_values = F[:S]
print(singular_values)
normCon = sum(singular_values)
normalise(x) = x/normCon
normalisedSVs = normalise.(singular_values)
print(normalisedSVs)
principal_components = F[:U] * diagm(F[:S])
#=
#figure(1, figsize=(4, 3))
clf()
#axes([.2, .2, .7, .7])
plot(singular_values, linewidth=2)
axis("tight")
xlabel("n components")
ylabel("explained variance")
=#

pyplot() # Switch to using the PyPlot.jl backend
#plotly() # alternatively, for backend for web interactivity
# Todo: PLOT DOES NOT MAKE SENSE, VALUES SHOULD BE MONOTONIC DECREASING
Plots.plot(singular_values, linewidth=2, xlabel = "principal components", ylabel = "Singular_values", title="Singular values")
Plots.bar!(singular_values, linewidth=1, xlim = [0, Inf], ylim=[0,2000], xlabel = "principal components", ylabel = "Singular_values", title="Singular values")
Plots.bar(normalisedSVs, linewidth=1, xlim = [0, Inf], ylim=[0,1], xlabel = "principal components", ylabel = "Explained variance", title="Singular values")
#=
# alternatively using Plots
using Plots
plotly() # Choose the Plotly.jl backend for web interactivity
plot(singular_values, linewidth=2, title="Singular values")
=#
pyplot() # Switch to using the PyPlot.jl backend
Plots.scatter(principal_components, linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components")
# .. shows that the first PC (blue) is very distinctive to the rest of the PCs
print(minimum(principal_components[1,:]))
print(minimum(principal_components[1,2:end])) # second smalles value in first Pc
print(maximum(principal_components[1,:]))
print(minimum(principal_components[2,:]))
print(maximum(principal_components[2,:]))

Plots.scatter(principal_components[1,:], principal_components[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
Plots.scatter(principal_components[1, 2:end], principal_components[2, 2:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
legend()
title("Projection onto 2 PCs")
savefig("Projection_2PCs")

# compute the centroid of the data points
x = principal_components[1,:]
y = principal_components[2,:]
centroid = (sum(x) / length(x), sum(y) / length(y))
print(centroid)

print(F[:S])
