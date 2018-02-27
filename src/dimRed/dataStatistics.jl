# -*- coding: utf-8 -*-
#=

The following program visualises the data set.
It computes and plots the mean and variance of expression levels of the 96 genes over the 547 cells;
it also computes the correlation matrix for two genes.

Note: You have to run Readin.jl first (dataframe data needed).

Args:
    -

Returns:
    -
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

# COMPUTES THE MEANS AND VARIANCES OF THE GENES IN DATA
function plotStatistics(data_array)
    means = []
    variance = []
    for col in 1:last(size(data_array))
        push!(means, mean(data_array[:, col]))
        push!(variance, var(data_array[:, col]))
    end
    #Plots.scatter(means)
    return means, variance
end

statistics = @time plotStatistics(data_array)
print("mean: ", statistics[1], "\n variance:", statistics[2])

# PLOTS MEAN AND VARIANCES OF THE GENES IN DATA
#Plots.plot(statistics[2], ylabel = "gene expression level", xlabel = "genes (1-96)", title = "mean/ variance of gene expression level", label = "mean", legend=true)
#Plots.plot!(statistics[1], label = "variance", legend = true)
Plots.bar(statistics[2], ylabel = "gene expression level", xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels", label = "mean", legend=true)
Plots.bar!(statistics[1], label = "variance", legend = true)
#Plots.scatter(statistics[1])
#Plots.scatter(statistics[2])
legend()
title("Mean / variance of data")
savefig("MeanVar_data")

meansF = Array{Float64}(statistics[1]);
varsF = Array{Float64}(statistics[2]);

names = DataFrames.names(data);
df_stats = DataFrame(gene_name = names, mean = statistics[1], variance = statistics[2])
print("Data statistics summary", df_stats)

# GENES SORTED BY DECREASING MEAN OF GENE EXPRESSION LEVEL
df_stats_Msort = sort(df_stats, cols = (:mean),
                    rev = (true));
# GENES SORTED by INCREASING VARIANCE, IE UNCERTAINTY
df_stats_Vsort = sort(df_stats, cols = cols = (:variance),
                    rev = (false));
# PLOTS
#@df data plot(1:548, [names], xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels, sorted by mean", legend=true)
@df df_stats_Msort plot(1:96, [:mean :variance], colour = [:blue :red], xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels, sorted by mean", legend=true)
@df df_stats_Vsort plot(1:96, [:mean :variance], colour = [:blue :red], xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels, sorted by variance", legend=true)
#@df df_stats_Msort bar(1:96, [:mean], colour = [:blue], xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels", legend=true)
#@df df_stats_Msort bar!(1:96, [:variance], colour = [:red], xlabel = "genes (1-96)", title = "mean/ variance of gene expression levels", legend=true)

#@df data corrplot(cols(1:96), grid = false)


###################
# TODO: BOXPLOTS
##################
#=
#@df df_stats_Msort violin(:mean,:variance,marker=(0.2,:blue,stroke(0)))
@df df_stats_Msort plot(violin(:mean, :value))
@df data boxplot(group=:names)
    @df data boxplot(:n)
=#
# maybe iterate over list of symbols
# THIS SHOULD BE WORKING!
#@df data boxplot(df_stats_Msort[1:3,1])
@df data boxplot(:Gapdh)
#@df data boxplot(:Actb, :Gapdh)
@df data boxplot!(:Actb)
@df data boxplot!(:Kdm1a, legend = true)



#@df df_stats_Msort violin(:mean, :gene_name)
#@df singers boxplot!(:VoicePart,:Height,marker=(0.3,:orange,stroke(2)))
#violin(df_stats[findin(df_stats[:gene_name],names[1:12]),:],:gene_names,:value, xlabel="", title="ViolinPlot of the gene expression data")

#= # from https://plot.ly/julia/box-plots/
using Plotly

x = (["day 1", "day 1", "day 1", "day 1", "day 1", "day 1",
      "day 2", "day 2", "day 2", "day 2", "day 2", "day 2"])


trace1 = [
  "y" => [0.2, 0.2, 0.6, 1.0, 0.5, 0.4, 0.2, 0.7, 0.9, 0.1, 0.5, 0.3],
  "x" => x,
  "name" => "kale",
  "marker" => ["color" => "#3D9970"],
  "type" => "box"]
trace2 = [
  "y" => [0.6, 0.7, 0.3, 0.6, 0.0, 0.5, 0.7, 0.9, 0.5, 0.8, 0.7, 0.2],
  "x" => x,
  "name" => "radishes",
  "marker" => ["color" => "#FF4136"],
  "type" => "box"]
trace3 = [
  "y" => [0.1, 0.3, 0.1, 0.9, 0.6, 0.6, 0.9, 1.0, 0.3, 0.6, 0.8, 0.5],
  "x" => x,
  "name" => "carrots",
  "marker" => ["color" => "#FF851B"],
  "type" => "box"]
data = [trace1, trace2, trace3]
layout = [
  "yaxis" => [
    "title" => "normalized moisture",
    "zeroline" => false
  ],
  "boxmode" => "group"]
response = Plotly.plot(data, ["layout" => layout, "filename" => "box-grouped", "fileopt" => "overwrite"])
plot_url = response["url"]
=#

###########
# GIVES HANDY DESCRIPTION / STATISTICS ON DATA FRAME
###########
#description = DataFrames.describe(data)
description = StatsBase.describe(data)
#description[:Wnt5a]

#names = names(data)

data[:Zfp281]

#=
# SELECT WHICH GENES YOU WANT TO HAVE PLOT
# GENERALISE
genes = DataFrames.names(data)
genes[1]
data[1]
gene1 = names[1]
gene2 = names[2]
@df data StatPlots.plot(1:1:547, :Actb, title="Expression of ")
#first(size(data_array)) # = 547
=#

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
