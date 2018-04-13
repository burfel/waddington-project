using StatPlots
using DataArrays
using Plots, GR
using PyPlot
gr()

sim = rand(5, 1000)

import PyPlot
x = 1:0.5:20
y = 1:0.5:10
f(x,y) = begin
        (3x + y ^ 2) * abs(sin(x) + cos(y))
    end
X = repmat(x',length(y),1)
Y = repmat(y,1,length(x))
Z = map(f,X,Y)
p1 = PyPlot.contour(x,y,f,fill=true)
p2 = PyPlot.contour(x,y,Z)
plot(p1,p2)


"""
# train a PCA model
M = fit(PCA, Xtr; maxoutdim=3)
P = projection(M)


indim(M)
principalvars(M)

typeof(P)
P(:,1)
scatter(P(:,1),P(:,2), P(3,:), title="My Scatter Plot")
"""
