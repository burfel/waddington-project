#=

Function that plots the landscape of a given data set.

Note:  - run Single_cell_data/reading_in_example.jl first
       - 3mostcorrelated.jl gives pairs of correlated genes which are plotted

Args:
    dataframe: data in the form of a dataframe or an array; here: 547 cells, 96 genes

Returns:
    -
=#


using Plots, DataStructures, DataFrames

function two_genes(dataframe)
           
           plotly()

           data_array = convert(Array,dataframe)

           data_array = collect(skipmissing(data_array))
           data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))
           data_array = data_array'

           # TAKES INPUT (2 genes that should be plottet)
           (no_genes,no_cells) = size(data_array)
           gene1 = input("Gene 1: ")
           gene2 = input("Gene 2: ")

           gene1 = convert(Symbol,gene1)
           gene2 = convert(Symbol,gene2)

           gene1index = find(names(data).==gene1)[1]
           gene2index = find(names(data).==gene2)[1]


           h=20
           Ps=zeros(eye(h))
           xspan=linspace(0.0001,20.0,h) #

           total = 0

           for i = 1:length(data_array[1,:])

                      sol = [data_array[gene1index,i],data_array[gene2index,i]]
                      x1=1
                      x2=1
                      found1=false
                      found2=false
                      impossible=false

                      while(!(found1 && found2))
                                 if (sol[1]<xspan[1] || sol[1]>=xspan[h] || sol[2]<xspan[1] || sol[2]>=xspan[h])
                                            impossible=true
                                            break;
                                 end

                                 if !(xspan[x1]<=sol[1] && sol[1]<xspan[x1+1]) && !found1
                                            x1=x1+1
                                 else
                                            found1=true
                                 end

                                 if !(xspan[x2]<=sol[2] && sol[2]<xspan[x2+1]) && !found2
                                            x2=x2+1
                                 else
                                            found2=true
                                 end
                      end

                      if !impossible
                                 Ps[x1,x2]=Ps[x1,x2]+1
                                 total=total+1

                      end
           end

           Ps=Ps/(total*(4/h)^2)+0.00001

           U=-log.(Ps)

           p1 = scatter(data[:,gene1index],data[:,gene2index],xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           p2 = surface(xspan,xspan,U,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])


           ####################################### Kernel Density Estimation method ############################################


           h=20
           Mtfinal=zeros(no_cells,2)
           xspan=linspace(0.5,20.0,h)

           total=0

           for i = 1:no_cells
               sol = [data_array[gene1index,i],data_array[gene2index,i]]
               Mtfinal[i,1]=sol[1]
               Mtfinal[i,2]=sol[2]

           end

           hps=1.06*1000^(-1/5)

           function fPs(x1,x2)

               sum = 0
               for i = 1:no_cells
                   sum=sum+exp(-((x1-Mtfinal[i,1])^2+(x2-Mtfinal[i,2])^2)/(2*(1.06*1000^(-1/5))^2))
               end

               sum=-log(sum/(1000*hps*sqrt(2*pi)))

               return sum
           end


           Results=zeros(h,h)

           for i = 1:h
               for j = 1:h
                   Results[i,j]=fPs(xspan[i],xspan[j])
               end
           end


           ########################################################################################################################


           p3 = surface(xspan,xspan,Results,xlabel=names(data)[gene1index],ylabel=names(data)[gene2index])

           return(p1,p2,p3)
end


# CALL THE FUNCTION

# ...on the whole data set
(p1,p2,p3) = two_genes(data) # Fgf4 and Gdf3
plot(p1)
plot(p2)
plot(p3)

# ...on t24 only
(p1,p2,p3) = two_genes(t24_only) # Pou5f1 and Trp53
plot(p1)
plot(p2)
plot(p3)

# ...on embryonic stem cells only
(p1,p2,p3) = two_genes(ESC_only) # Dnmt3b and Prmt7
plot(p1)
plot(p2)
plot(p3)

#most correlated - good
(Fgf4_Gdf3_p1,Fgf4_Gdf3_p2,Fgf4_Gdf3_p3) = two_genes(data)
#most entropy - rubbish
(Actb_Gapdh_p1,Actb_Gapdh_p2,Actb_Gapdh_p3) = two_genes(data)
#most mutual information - good
(Fgf4_Pou5f1_p1,Fgf4_Pou5f1_p2,Fgf4_Pou5f1_p3) = two_genes(data)

#Nanog and Gata6 were plotted in RDB presentation
(Nanog_Gata6_p1,Nanog_Gata6_p2,Nanog_Gata6_p3) = two_genes(data)
plot(Nanog_Gata6_p3)

(x,y,z1) = two_genes(t24_only)
plot(z1, title = "24 hours")

(x,y,z2) = two_genes(t48_only)
plot(z2, title = "48 hours")

(x,y,z3) = two_genes(t72_only)
plot(z3, title = "72 hours")

(x,y,z4) = two_genes(t168_only)
plot(z4,title = "168 hours")
