#= DEPRECATED 
The following program plots the explained variance of the data points in data_array
Note: You have to run topCorr.jl first (dataframe data needed).
=#

using Plots
using PyPlot

# SVD without ScikitLearn:
data_arrayT = transpose(data_array)
FT = svdfact(data_arrayT)
#F[:U]
#F[:S]
#F[:V]
singular_valuesT = FT[:S]
print(singular_valuesT)
normConT = sum(singular_valuesT)
normalise(x) = x/normConT
normalisedSVsT = normalise.(singular_valuesT)
print(normalisedSVsT)
principal_componentsT = FT[:U] * diagm(FT[:S])
#principal_components = diagm(F[:S]) * F[:Vt]
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
Plots.plot(singular_valuesT, linewidth=2, xlabel = "principal components", ylabel = "Singular_values", title="Singular values")
Plots.bar!(singular_valuesT, linewidth=1, xlim = [0, Inf], ylim=[0,2000], xlabel = "principal components", ylabel = "Singular_values", title="Singular values")
#legend()
#title("Singular values")
#savefig("singular_values_of_resp_PCs")
Plots.bar(normalisedSVsT, linewidth=1, xlim = [0, Inf], ylim=[0,1], xlabel = "principal components", ylabel = "Explained variance", title="Singular values")
#legend()
#title("Explained variance of PCs")
#savefig("singular_values_of_resp_PCs_normed")

# Todo: maybe cumulative distrib of explained variance here

#=
# alternatively using Plots
using Plots
plotly() # Choose the Plotly.jl backend for web interactivity
plot(singular_values, linewidth=2, title="Singular values")
=#
pyplot() # Switch to using the PyPlot.jl backend
Plots.scatter(principal_componentsT, linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components")
# .. shows that the first PC (blue) is very distinctive to the rest of the PCs
print(minimum(principal_componentsT[1,:]))
print(minimum(principal_componentsT[1,2:end])) # second smalles value in first Pc
print(maximum(principal_componentsT[1,:]))
print(minimum(principal_componentsT[2,:]))
print(maximum(principal_componentsT[2,:]))

# to get transformed data, we multiply by V; this equals the principal components U*S
to_new_basisT = data_arrayT * FT[:V]


Plots.scatter(to_new_basisT[1,:], to_new_basisT[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
# equals...
Plots.scatter(principal_componentsT[1,:], principal_componentsT[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
#Plots.scatter(principal_components[1, 2:end], principal_components[2, 2:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")
#legend()
#title("Projection onto 2 PCs")
#savefig("Projection_2PCs")

# compute the centroid of the data points
x = principal_componentsT[1,:]
y = principal_componentsT[2,:]
centroid = (sum(x) / length(x), sum(y) / length(y))
print(centroid)

"""
outputT = data_array * principal_componentsT
print(outputT)
Plots.scatter(output)

#principal_components[1,1].*data_array[:,1]
"""
