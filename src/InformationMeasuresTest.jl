using InformationMeasures

function df2array(dataframe)
       data_array = convert(Array,dataframe)
       data_array = collect(skipmissing(data_array))
       data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))
       return data_array
end

size(df2array(data))

function maxentropy(dataframe)
       data_array = df2array(dataframe)

       (cells,genes) = size(data_array)
       ent = zeros(genes,genes)
       entmax = 0.0
       index1 = 1
       index2 = 1
       for i = 1:genes
              for j = i:genes
                     ent[i,j] = get_entropy(data_array[:,i],data_array[:,j])
                     ent[j,i] = ent[i,j]
                     if ent[i,j]>entmax
                            entmax = ent[i,j]
                            (index1, index2) = (i,j)
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2]
end

maxentropy(data)
maxentropy(t24_only)
maxentropy(ESC_only)

function maxconditionalentropy(dataframe)
       data_array = df2array(dataframe)

       (cells,genes) = size(data_array)
       ent = zeros(genes,genes)
       entmax = 0.0
       index1 = 1
       index2 = 1
       for i = 1:genes
              for j = i:genes
                     ent[i,j] = get_conditional_entropy(data_array[:,i],data_array[:,j])
                     ent[j,i] = ent[i,j]
                     if ent[i,j]>entmax
                            entmax = ent[i,j]
                            (index1, index2) = (i,j)
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2]
end

maxconditionalentropy(data)
maxconditionalentropy(t24_only)
maxconditionalentropy(ESC_only)

function maxmutualinformation(dataframe)
       data_array = df2array(dataframe)

       (cells,genes) = size(data_array)
       ent = zeros(genes,genes)
       entmax = 0.0
       index1 = 1
       index2 = 1
       for i = 1:genes
              for j = i+1:genes
                     ent[i,j] = get_mutual_information(data_array[:,i],data_array[:,j])
                     ent[j,i] = ent[i,j]
                     if ent[i,j]>entmax
                            entmax = ent[i,j]
                            (index1, index2) = (i,j)
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2]
end

maxmutualinformation(data)
maxmutualinformation(t24_only)
maxmutualinformation(ESC_only)

function maxinteractioninformation(dataframe)
       data_array = df2array(dataframe)

       (cells,genes) = size(data_array)
       ent = zeros(genes,genes,genes)
       entmax = 0.0
       index1 = 1
       index2 = 1
       index3 = 1
       for i = 1:genes
              for j = 1:genes
                     for k = 1:genes
                            ent[i,j,k] = get_interaction_information(data_array[:,i],data_array[:,j],data_array[:,k])
                            if ent[i,j,k]>entmax
                                   entmax = ent[i,j,k]
                                   (index1, index2, index3) = (i,j,k)
                            end
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2], names(dataframe)[index3]
end

maxinteractioninformation(data)
