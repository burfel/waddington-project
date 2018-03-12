# -*- coding: utf-8 -*-
#=

Function that computes a generalised PCA (SVD).
It plots the explained variance of the data points in the input data set.
Note: You have to run Readin.jl first (dataframe data needed).

Args:
    data_array: data in form of an array, here: 547 x 96 matrix

Returns:
    The principal component matrix.

=#


#using ScikitLearn
# NOT GOOD: LOADING MULTIPLE PLOTTING PACKAGES
using Plots
using PyPlot
using StatPlots
using DataFrames
using StatsBase
#using ScikitLearn.GridSearch: GridSearchCV
#using ScikitLearn.Pipelines: Pipeline, named_steps

#@sk_import decomposition: PCA
#@sk_import datasets: load_digits
#@sk_import linear_model: LogisticRegression


function generalPCA(data_array)

# SVD without ScikitLearn:

    # Convert the data if required
    if data_array isa DataFrame
        data_array = convert(Array,data_array)
    elseif data_array isa Array
        data_array = data_array
    else
        print("Data neither of type DataFrame or Array")
    end

    F = svdfact(data_array)
    #F[:U]
    #F[:S]
    #F[:V]
    singular_values = F[:S]
    print("Singular values: ", singular_values)

    normCon = sum(singular_values)
    normalise(x) = x/normCon
    normalisedSVs = normalise.(singular_values)
    print("Normalised singular values: ", normalisedSVs)
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

    #  PLOTS
    pyplot() # Switch to using the PyPlot.jl backend
    #plotly() # alternatively, for backend for web interactivity

    Plots.plot(singular_values, linewidth=2, xlabel = "principal components", ylabel = "Singular values", title="Singular values", legend=false)
    Plots.bar!(singular_values, linewidth=1, xlim = [0, Inf], ylim=[0,2000], xlabel = "principal components", ylabel = "Singular values", title="Singular values", label=false)
    #PyPlot.legend()
    #PyPlot.title("Singular values")
    #PyPlot.savefig("singular_values_of_resp_PCs")

    Plots.plot(normalisedSVs, linewidth=2, xlabel = "principal components", ylabel = "Singular values", title="Singular values", legend=false)
    Plots.bar!(normalisedSVs, linewidth=1, xlim = [0, Inf], ylim=[0,1], xlabel = "principal components", ylabel = "Explained variance", title="Singular values", legend=false)
    #PyPlot.legend()
    #PyPlot.title("Explained variance of PCs")
    #PyPlot.savefig("singular_values_of_resp_PCs_normed")


    # TODO: maybe cumulative distrib of explained variance here

    #=
    # alternatively using Plots
    using Plots
    plotly() # Choose the Plotly.jl backend for web interactivity
    plot(singular_values, linewidth=2, title="Singular values")
    =#
    pyplot() # Switch to using the PyPlot.jl backend
    # not very meaningful plot: # shows that first PC is a linear combination of many many cells, rest can be explained by single cells
    Plots.scatter(transpose(principal_components), linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components", legend=false)
    Plots.scatter(principal_components, linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components", legend=false)
    # .. shows that the first PC (blue) is very distinctive to the rest of the PCs
    #print(minimum(principal_components[1,:]))
    #print(minimum(principal_components[1,2:end])) # second smalles value in first Pc
    #print(maximum(principal_components[1,:]))
    #print(minimum(principal_components[2,:]))
    #print(maximum(principal_components[2,:]))

    # cells projected onto two first components
    Plots.scatter(principal_components[1,:], principal_components[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
    Plots.scatter(principal_components[1, 2:end], principal_components[2, 2:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
    #PyPlots.legend()
    #PyPlots.title("Projection onto 2 PCs")
    #PyPlots.savefig("Projection_2PCs")

    # COMPUTE THE CENTROID OF THE DATA POINTS
    x = principal_components[1,:]
    y = principal_components[2,:]
    centroid = (sum(x) / length(x), sum(y) / length(y))
    print("Centroid: ", centroid)

    # to get transformed data, we multiply by V; this equals the principal components U*S
    output = principal_components

    return output

end

generalPCA(data_array)
