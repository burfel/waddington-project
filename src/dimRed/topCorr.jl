# -*- coding: utf-8 -*-
using DataFrames

"""Function that selects the most correlated genes.
Note: -

Args:
    data: data (cells | genes) as a csv file; here: 547 cells, 96 genes
    amount (optional): number of top genes one wants to get, default: 3

Returns:
    The top #amount correlated genes.
"""

df1= readtable("../Single_cell_data/ND-data.txt",separator=' ',header=true);
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
data = delete!(df1, :variable)

ann=readtable("../Single_cell_data/annotation.txt",header=true,separator=' ');

ESC_only = data[ann[:ThreeST] .== "ESC", :]

# only data sampled at 24 hours
t24_only = data[ann[:Time] .== 24, :]
