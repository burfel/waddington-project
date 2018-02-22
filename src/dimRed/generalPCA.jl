#=
The following program computes a generalised PCA (SVD).
It plots the explained variance of the data points in the input data set.

Note: You have to run topCorr.jl first (dataframe data needed).
=#

#using ScikitLearn
# NOT GOOD: LOADING MULTIPLE PLOTTING PACKAGES
using Plots
using PyPlot
using StatPlots
using DataFrames
#using ScikitLearn.GridSearch: GridSearchCV
#using ScikitLearn.Pipelines: Pipeline, named_steps

#@sk_import decomposition: PCA
#@sk_import datasets: load_digits
#@sk_import linear_model: LogisticRegression

# COMPUTES THE MEANS AND VARIANCES OF THE GENES IN DATA
function plotStatistics(data_array)
    means = []
    variance = []
    for col in 1:last(size(data_array))
        push!(means, mean(data_array[:, col]))
        push!(variance, var(data_array[:, col]))
    end
    Plots.scatter(means)
    return means, variance
end

@time plotStatistics(data_array)

statistics = plotStatistics(data_array)
print("mean: ", statistics[1], "\n variance:", statistics[2])
# Todo: PRINT THAT NICER IN A DATAFRAME


# PLOTS MEAN AND VARIANCES OF THE GENES IN DATA
Plots.plot(statistics[2], ylabel = "gene expression level", xlabel = "genes (1-96)", title = "mean/ variance of gene expression level", label = "mean", legend=true)
Plots.plot!(statistics[1], label = "variance", legend = true)
#Plots.scatter(statistics)

meansF = Array{Float64}(statistics[1])
varsF = Array{Float64}(statistics[2])
# not working:
#convert(DataFrame, meansF)
# Todo: convert to dataframe, sort, plot

#showcols(data)

# GIVES HANDY DESCRIPTION / STATISTICS ON DATA FRAME
description = DataFrames.describe(data)
#description[:Wnt5a]

#names = names(data)

data[:Zfp281]

# SELECT WHICH GENES YOU WANT TO HAVE PLOT
# GENERALISE
genes = DataFrames.names(data)
genes[1]
data[1]
gene1 = names[1]
gene2 = names[2]
@df data StatPlots.plot(1:1:547, :Actb, title="Expression of ")
#first(size(data_array)) # = 547

#=
println("Column\tMeanX\tMedianX\tStdDev X\tMeanY\t\t\tStdDev Y\t\tCorr\t")
map(xcol -> println(
        xcol,                   "\t",
        mean(anscombe[xcol]),   "\t",
        median(anscombe[xcol]), "\t",
        std(anscombe[xcol]),    "\t",
        mean(anscombe[xcol]),   "\t",
        std(anscombe[xcol]),    "\t",
        cor(anscombe[xcol], anscombe[xcol])),

    [:X1, :X2, :X3, :X4],
    [:Y1, :Y2, :Y3, :Y4]);
=#
#=
println("Column\tMean\tMedian\tVariance\tStdDev\tCovariance\tCorrelation\t")

map(xcol -> println(
    xcol,               "\t",
    mean(data[xcol]),   "\t",
    median(data[xcol]), "\t",
    var(data[xcol]),    "\t",
    std(data[xcol]),    "\t",
    cov(data[xcol]),    "\t",
    cor(data[xcol]),

    [names];
=#


function generalPCA(data_array)

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
    legend()
    title("Singular values")
    savefig("singular_values_of_resp_PCs")
    Plots.bar(normalisedSVs, linewidth=1, xlim = [0, Inf], ylim=[0,1], xlabel = "principal components", ylabel = "Explained variance", title="Singular values")
    legend()
    title("Explained variance of PCs")
    savefig("singular_values_of_resp_PCs_normed")

    # Todo: maybe cumulative distrib of explained variance here

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

    # to get transformed data, we multiply by V; this equals the principal components U*S
    output = principal_components

    return output

end
