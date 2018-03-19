# -*- coding: utf-8 -*-
#=

Function that plots the landscape of two chosen genes from a given data set.

Note:  - run Single_cell_data/reading_in_example.jl first
       - run simulations/PFM_reduction_2D.jl first
       - run simulations/KDE_reduction_2D.jl first
       - 3mostcorrelated.jl gives pairs of correlated genes which are plotted

Args:
    dataframe: data in the form of a dataframe or an array; here: 547 cells, 96 genes

Returns:
    -p1 a scatter of the data
    -p2 a plot of the landscape given by the Probability Flux method
    -p3 a plot of the landscape given by the Kernel Density Estimation method
=#


using Plots, DataStructures, DataFrames
plotlyjs()
function two_genes(dataframe)
           data_array = convert(Array,dataframe)
           data_array = collect(skipmissing(data_array))
           data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))

           # TAKES INPUT (2 genes that should be plotted)
           (no_genes,no_cells) = size(data_array)
           gene1 = input("Gene 1: ")
           gene2 = input("Gene 2: ")

           gene1 = convert(Symbol,gene1)
           gene2 = convert(Symbol,gene2)

           gene1index = find(names(data).==gene1)[1]
           gene2index = find(names(data).==gene2)[1]

           s = scatter(data[:,gene1index],data[:,gene2index],xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           ###

           end_state = [data_array[:,gene1index] data_array[:,gene2index]]
           boundaries = [0.0001 0.0001; 20.0 20.0]
           nb_bin = 10
           dim1 = 1
           dim2 = 2

           xspan = Span_2D(boundaries,nb_bin,dim1,dim2)

           #U = PFM_reduction_2D(end_state,boundaries,nb_bin,dim1,dim2)
           #p2 = surface(xspan[1,:],xspan[2,:],U,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           U_kde = KDE_reduction_2D(end_state, boundaries, nb_bin, dim1,dim2)
           U = surface(xspan[1,:],xspan[2,:],U_kde,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           return(s,U)
end

function two_genes_terminal(dataframe)
           data_array = convert(Array,dataframe)
           data_array = collect(skipmissing(data_array))
           data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))

           # TAKES INPUT (2 genes that should be plotted)
           (no_genes,no_cells) = size(data_array)
           print("Gene 1: ")
           gene1 = readline(STDIN)
           print("Gene 2: ")
           gene2 = readline(STDIN)

           gene1 = convert(Symbol,gene1)
           gene2 = convert(Symbol,gene2)

           gene1index = find(names(data).==gene1)[1]
           gene2index = find(names(data).==gene2)[1]

           s = scatter(data[:,gene1index],data[:,gene2index],xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           ###

           end_state = [data_array[:,gene1index] data_array[:,gene2index]]
           boundaries = [0.0001 0.0001; 20.0 20.0]
           nb_bin = 20
           dim1 = 1
           dim2 = 2

           xspan = Span_2D(boundaries,nb_bin,dim1,dim2)

           #U = PFM_reduction_2D(end_state,boundaries,nb_bin,dim1,dim2)
           #p2 = surface(xspan[1,:],xspan[2,:],U,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           U_kde = KDE_reduction_2D(end_state, boundaries, nb_bin, dim1,dim2)
           U = surface(xspan[1,:],xspan[2,:],U_kde,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           return(s,U)
end

function two_genes_input(dataframe,gene1,gene2)

           data_array = convert(Array,dataframe)
           data_array = collect(skipmissing(data_array))
           data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))

           # TAKES INPUT (2 genes that should be plotted)
           (no_genes,no_cells) = size(data_array)
           #gene1 = input("Gene 1: ")
           #gene2 = input("Gene 2: ")

           gene1 = convert(Symbol,gene1)
           gene2 = convert(Symbol,gene2)

           gene1index = find(names(data).==gene1)[1]
           gene2index = find(names(data).==gene2)[1]

           s = scatter(data[:,gene1index],data[:,gene2index],xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           ###

           end_state = [data_array[:,gene1index] data_array[:,gene2index]]

           boundaries = [0.0001 0.0001; 20.0 20.0]
           nb_bin = 20
           dim1 = 1
           dim2 = 2

           xspan = Span_2D(boundaries,nb_bin,dim1,dim2)

           #U = PFM_reduction_2D(end_state,boundaries,nb_bin,dim1,dim2)
           #p2 = surface(xspan[1,:],xspan[2,:],U,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           U_kde = KDE_reduction_2D(end_state, boundaries, nb_bin, dim1,dim2)
           U = surface(xspan[1,:],xspan[2,:],U_kde,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           return(s,U)
end

function two_genes_skip_zeros(dataframe,gene1,gene2)

           data_array = convert(Array,dataframe)
           data_array = collect(skipmissing(data_array))
           data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))

           # TAKES INPUT (2 genes that should be plotted)
           (no_genes,no_cells) = size(data_array)
           #gene1 = input("Gene 1: ")
           #gene2 = input("Gene 2: ")

           gene1 = convert(Symbol,gene1)
           gene2 = convert(Symbol,gene2)

           gene1index = find(names(data).==gene1)[1]
           gene2index = find(names(data).==gene2)[1]

           s = scatter(data[:,gene1index],data[:,gene2index],xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           ###

           end_state = [data_array[:,gene1index] data_array[:,gene2index]]

           ### SKIP ZEROS
           (rows,cols)=size(end_state)
           nzs=[]
           count=1
           for n=1:rows
                      (i,j) = end_state[n,:]
                      if !(i==0.00 || j==0.00)
                                 push!(nzs,end_state[n,:])
                                 count+=1
                      end
           end
           nonzeros= zeros(length(nzs),2)
           count=1
           for n=1:rows
                      (i,j) = end_state[n,:]
                      if !(i==0.00 || j==0.00)
                                 nonzeros[count,:]=end_state[n,:]
                                 count+=1
                      end
           end

           boundaries = [0.0001 0.0001; 20.0 20.0]
           nb_bin = 20
           dim1 = 1
           dim2 = 2

           xspan = Span_2D(boundaries,nb_bin,dim1,dim2)

           #U = PFM_reduction_2D(nonzeros,boundaries,nb_bin,dim1,dim2)
           #p2 = surface(xspan[1,:],xspan[2,:],U,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           U_kde = KDE_reduction_2D(nonzeros, boundaries, nb_bin, dim1,dim2)
           U = surface(xspan[1,:],xspan[2,:],U_kde,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           return(s,U)
end


## CHOOSE BACKEND

## CALL THE FUNCTION

## ...on the whole data set
# (p1,p2,p3) = two_genes(data) # Fgf4 and Gdf3
# plot(p1)
# plot(p2)
# plot(p3)
#
## ...on t24 only
# (p1,p2,p3) = two_genes(t24_only) # Pou5f1 and Trp53
# plot(p1)
# plot(p2)
# plot(p3)
#
## ...on embryonic stem cells only
# (p1,p2,p3) = two_genes(ESC_only) # Dnmt3b and Prmt7
# plot(p1)
# plot(p2)
# plot(p3)
#
# #most correlated - good
# (Fgf4_Gdf3_p1,Fgf4_Gdf3_p2,Fgf4_Gdf3_p3) = two_genes(data)
# #most entropy - rubbish
# (Actb_Gapdh_p1,Actb_Gapdh_p2,Actb_Gapdh_p3) = two_genes(data)
# #most mutual information - good
# (Fgf4_Pou5f1_p1,Fgf4_Pou5f1_p2,Fgf4_Pou5f1_p3) = two_genes(data)
#
# #Nanog and Gata6 were plotted in RDB presentation
# (Nanog_Gata6_p1,Nanog_Gata6_p2,Nanog_Gata6_p3) = two_genes(data)
# plot(Nanog_Gata6_p3)
#
# (x,y,z1) = two_genes(t24_only)
# plot(z1, title = "24 hours")
#
# (x,y,z2) = two_genes(t48_only)
# plot(z2, title = "48 hours")
#
# (x,y,z3) = two_genes(t72_only)
# plot(z3, title = "72 hours")
# #
# (x,y,z4) = two_genes(t168_only)
#  plot(z4,title = "168 hours")
#  #
