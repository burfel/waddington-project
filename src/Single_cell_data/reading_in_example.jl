# -*- coding: utf-8 -*-
#=

Function that reads in the data set ND-data.txt and stores it in a dataframe.
Note: You can modify the program  so that only selected data is read in;
      annotation.txt may also be modified to provide extra information.

=#

using DataFrames

# Reading in the data from the ND-data text file to a dataframe
df1= readtable("ND-data.txt",separator=' ',header=true);
cellnames=["Names"]
for i in 1:547
    if i<10
    push!(cellnames,"cell00$i")
    elseif i<100
        push!(cellnames,"cell0$i")
    else
        push!(cellnames,"cell$i")
    end
end
names!(df1,[Symbol(X) for X in cellnames]);
df1=stack(df1,names(df1)[2:end]);
df1=unstack(df1,:variable,:Names,:value)

# Delete the column containing the cell name.
# An extra step from my classification work before, but perhaps still helpful
data = delete!(df1, :variable)

# Extra information provided by the annotation.txt file
ann=readtable("annotation.txt",header=true,separator=' ');


# Then you can do things like get only ESC data
ESC_only = data[ann[:ThreeST] .== "ESC", :]

EPI_only = data[ann[:ThreeST] .== "EPI", :]

NPC_only = data[ann[:ThreeST] .== "NPC", :]
# ... or only data sampled at 24 hours
t24_only = data[ann[:Time] .== 24, :]

t48_only = data[ann[:Time] .== 48, :]

t72_only = data[ann[:Time] .== 72, :]

t96_only = data[ann[:Time] .== 96, :]

t120_only = data[ann[:Time] .== 120, :]

t144_only = data[ann[:Time] .== 144, :]

t168_only = data[ann[:Time] .== 168, :]

t0_only = data[ann[:Time] .== 0, :]
