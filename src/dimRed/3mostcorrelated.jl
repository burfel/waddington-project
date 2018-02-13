# 547 cells, 96 genes
using DataFrames

function top3correlated(data)
    # Convert the data if required
    if data isa DataFrame
        data = convert(Array,data)
    elseif data isa Array

    else
        print("Data neither of type DataFrame or Array")

    end

    (cells,genes)=size(data)

    CorMatrix = cor(data)
    for i=1:genes
        CorMatrix[i,i]=0.0
    end


    entry=0.0
    most_correlated1 = [0,0]
    for i=1:genes
        for j=i+1:genes
            if abs(CorMatrix[i,j]) > entry
                entry = CorMatrix[i,j]
                most_correlated1[1] = i
                most_correlated1[2] = j
            else
            end
        end
    end

    print(most_correlated1)

    CorMatrix[most_correlated1[1],most_correlated1[2]] = 0.0

    entry=0.0
    most_correlated2 = [0,0]
    for i=1:genes
        for j=i+1:genes
            if abs(CorMatrix[i,j]) > entry
                entry = CorMatrix[i,j]
                most_correlated2[1] = i
                most_correlated2[2] = j
            else
            end
        end
    end

    print(most_correlated2)




end


top3correlated(data)
top3correlated(t24_only)
top3correlated(ESC_only) # embryonic stem cells only


using Plots
scatter(data[:,22],data[:,31],data[:,67],xlabel=names(data)[22],ylabel=names(data)[31],zlabel=names(data)[67])
scatter(t24_only[:,67],t24_only[:,88],t24_only[:,90],xlabel=names(data)[67],ylabel=names(data)[88],zlabel=names(data)[90])
