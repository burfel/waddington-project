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
