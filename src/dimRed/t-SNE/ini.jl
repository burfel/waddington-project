Pkg.clone("https://github.com/zhmz90/BHTsne.jl.git")# installation

using BHTsne
using RDatasets
using Gadfly

iris = dataset("datasets", "iris")
samples = convert(Array, iris[:,1:4])
labels  = convert(Array, iris[:,5])

results = bh_tsne(samples, perplexity=30, verbose=true)

p = plot(x=results[:,1], y=results[:,2], color=labels)
draw(PNG("tsne_of_iris.png", 8inch, 6inch), p)

function bh_tsne(samples; no_dims=2, initial_dims=50, perplexity=50, theta=0.5, randseed=-1, verbose=false)
