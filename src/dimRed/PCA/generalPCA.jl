# as normal code, not in a function

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
#PyPlots.legend()
PyPlot.title("Singular values")
PyPlot.savefig("singular_values_of_resp_PCs")

Plots.plot(normalisedSVs, linewidth=2, xlabel = "principal components", ylabel = "Singular values", title="Singular values", legend=false)
Plots.bar!(normalisedSVs, linewidth=1, xlim = [0, Inf], ylim=[0,1], xlabel = "principal components", ylabel = "Explained variance", title="Singular values", legend=false)
#PyPlots.legend()
PyPlot.title("Explained variance of PCs")
PyPlot.savefig("singular_values_of_resp_PCs_normed")#PyPlots.legend()


# TODO: maybe cumulative distrib of explained variance here

#=
# alternatively using Plots
using Plots
plotly() # Choose the Plotly.jl backend for web interactivity
plot(singular_values, linewidth=2, title="Singular values")
=#
pyplot() # Switch to using the PyPlot.jl backend
# plot not very meaningful

Plots.scatter(principal_components, linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components", legend=false)
# .. shows that the first PC (blue) is very distinctive to the rest of the PCs
#print(minimum(principal_components[1,:]))
#print(minimum(principal_components[1,2:end])) # second smalles value in first Pc
#print(maximum(principal_components[1,:]))
#print(minimum(principal_components[2,:]))
#print(maximum(principal_components[2,:]))

######--------

# Visualise position of first two PCs
Plots.scatter(principal_components[1,:], principal_components[2,:], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Position of the 96 PCs")
Plots.scatter(principal_components[1, 2:end], principal_components[2, 2:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Position of the last 95 PCs")
Plots.scatter(principal_components[1, 3:end], principal_components[2, 3:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Position of the last 94 PCs")
Plots.scatter(principal_components[1, 4:end], principal_components[2, 4:end], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Position of the last 93 PCs")

Plots.scatter(principal_components[1,:], linewidth=2, xlabel = "Principal components", ylabel = "cells", label = "cell 1", title="")
Plots.scatter!(principal_components[2,:], linewidth=2, xlabel = "Principal components", ylabel = "cells", label = "cell 2", title="")
Plots.scatter!(principal_components[3,:], linewidth=2, xlabel = "Principal components", ylabel = "cells", label = "cell 3", title="")

# or all:
# not very meaningful plot: # shows that first PC is a linear combination of many many cells, rest can be explained by single cells
Plots.scatter(transpose(principal_components), linewidth=2, xlabel = "cells", ylabel ="gene expression of particular gene combination", title="Principal components", legend=false)


# Project data onto first two PCs
Plots.scatter(principal_components[:,1], principal_components[:,2], linewidth=2, xlabel = "PC1", ylabel = "PC2", title="Projection onto 2 PCs")

#PyPlots.legend()
#PyPlots.title("Projection onto 2 PCs")
#PyPlots.savefig("Projection_2PCs")

######-------

# COMPUTE THE CENTROID OF THE DATA POINTS
x = principal_components[1,:]
y = principal_components[2,:]
centroid = (sum(x) / length(x), sum(y) / length(y))
print("Centroid: ", centroid)
