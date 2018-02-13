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

# ... or only data sampled at 24 hours
t24_only = data[ann[:Time] .== 24, :]
