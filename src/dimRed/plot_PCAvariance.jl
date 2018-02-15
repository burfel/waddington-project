#=
The following program plots the explained variance of the data points in data_array
Note: You have to run topCorr.jl first (dataframe data needed).
=#

using ScikitLearn
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
# Todo: PLOT DOES NOT MAKE SENSE, VALUES SHOULD BE MONOTONIC DECREASING
histogram(singular_values, linewidth=2, xlabel = "principal components", ylabel = "Singular_values", title="Singular values")
histogram(normalisedSVs, linewidth=2, xlabel = "principal components", ylabel = "Explained variance", title="Singular values")
#=
# alternatively using Plots
using Plots
plotly() # Choose the Plotly.jl backend for web interactivity
plot(singular_values, linewidth=2, title="Singular values")
=#
pyplot() # Switch to using the PyPlot.jl backend
scatter(principal_components, linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components")
# .. shows that the first PC (blue) is very distinctive to the rest of the PCs
scatter(principal_components[1,:], principal_components[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
