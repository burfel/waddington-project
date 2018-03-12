# -*- coding: utf-8 -*-
#=

Function that reads in the data set given as a csv file and returns it as a dataframe or as an array.
Note: -

Args:
    loc: location (directory) of data (cells | genes) as a csv file; here: 547 cells, 96 genes
    df (optional): whether one wants it returned as a dataframe; default: array

Returns:
    The data as a dataframe or data array, depending on df.

=#

using DataFrames

function readin(loc, df=false)

#loc = ARGS[1]
    df1= readtable(loc*"ND-data.txt", separator=' ',header=true);
    #df1= readtable("../Single_cell_data/ND-data.txt",separator=' ',header=true);
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
    df1=DataFrames.stack(df1,DataFrames.names(df1)[2:end]);
    df1=DataFrames.unstack(df1,:variable,:Names,:value)

    # Delete the column containing the cell name.
    data = delete!(df1, :variable)

    #ann=readtable("../Single_cell_data/annotation.txt",header=true,separator=' ');
    ann=readtable(loc*"annotation.txt",header=true,separator=' ');

    ESC_only = data[ann[:ThreeST] .== "ESC", :]

    # only data sampled at 24 hours
    t24_only = data[ann[:Time] .== 24, :]


    data_array = convert(Array, data)
    #data_array = collect(skipmissing(data_array))
    data_array = reshape(data_array, (547,96))

    if df
        return data
    else
        return data_array
    end
end

#=
using Readin
Readin.readin("../Single_cell_data/", true)
=#

data = readin("../Single_cell_data/", true)
