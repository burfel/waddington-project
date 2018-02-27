# -*- coding: utf-8 -*-
#=

Program that performs different entropy, mutual information (MI) and higher order
measures computations from information theory on given data set.

1. Entropy
ent_1 = get_entropy(data_1)
ent_12 = get_entropy(data_1, data_2)
ent_123 = get_entropy(data_1, data_2, data_3)

2. Conditional entropy
ce_1_on_2 = get_conditional_entropy(data_1, data_2)

3. Mutual information
mi_12 = get_mutual_information(data_1, data_2)

4. Conditional mutual information
cmi_12_on_3 = get_conditional_mutual_information(data_1, data_2, data_3)

5. Interaction information
ii_123 = get_interaction_information(data_1, data_2, data_3)

6. Highest covariance

7. Total correlation
tc_123 = get_total_correlation(data_1, data_2, data_3)

8. Partial information decomposition
pid_123 = get_partial_information_decomposition(data_1, data_2, data_3)


Note:
    -

Args:
    -

Returns:
    -
=#

using InformationMeasures

# converts dataframe to array
function df2array(dataframe)
       data_array = convert(Array,dataframe)
       #data_array = collect(skipmissing(data_array))
       data_array = reshape(data_array,(nrow(dataframe),ncol(dataframe)))
       return data_array
end

size(df2array(data))

# 1. returns MAXIMAL (pairwise) ENTROPY
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

maxentrop = @time maxentropy(data)
maxentropy(t24_only)
maxentropy(ESC_only)

print("maximal entropy: ", maxentrop)


# 2. returns MAXIMAL CONDITIONAL ENTROPY
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

maxcondentrop = @time maxconditionalentropy(data)
maxconditionalentropy(t24_only)
maxconditionalentropy(ESC_only)

print("maximal conditional entropy: ", maxcondentrop)


# 3. returns MAXIMAL MUTUAL INFORMATION (MI) of two genes
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

maxMI = @time maxmutualinformation(data)
maxmutualinformation(t24_only)
maxmutualinformation(ESC_only)

print("maximal MI: ", maxMI)


# 4. returns CONDITIONAL MUTUAL INFORMATION of two genes --------------DOES NOT TERMINATE!
function maxCondMI(dataframe)
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
                            ent[i,j,k] = get_conditional_mutual_information(data_array[:,i],data_array[:,j],data_array[:,k])
                            if ent[i,j,k]>entmax
                                   entmax = ent[i,j,k]
                                   (index1, index2, index3) = (i,j,k)
                            end
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2], names(dataframe)[index3]
end

maxCondMInf = @time maxCondMI(data)

print("maximal conditional mutual information: ", maxCondMInf)


# 5. returns MAXIMAL INTERACTION INFORMATION of two genes --------------DOES NOT TERMINATE!
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

maxInteractionInfo = @time maxinteractioninformation(data)

print("maximal interaction information: ", maxInteractionInfo)


# 6. returns genes with HIGHEST VARIANCE
# TODO: maybe renaming function since name misleading....
function mostvaried(dataframe)
       data_array = df2array(dataframe)

       (cells,genes) = size(data_array)
       CovMat = zeros(genes,genes)
       maxcov = 0.0
       index1 = 1
       index2 = 1
       for i = 1:genes
              for j = i+1:genes
                     CovMat[i,j] = cov(data_array[:,i],data_array[:,j])
                     if CovMat[i,j]>maxcov
                            maxcov = CovMat[i,j]
                            index1 = i
                            index2 = j
                     end
              end
       end
       return maxcov, names(dataframe)[index1], names(dataframe)[index2]
end

highCov = @time mostvaried(data)

print("highest covariance: ", highCov)



# 7. returns TOTAL CORRELATION --------------DOES NOT TERMINATE!
function totalCorrelation(dataframe)
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
                            ent[i,j,k] = get_total_correlation(data_array[:,i],data_array[:,j],data_array[:,k])
                            if ent[i,j,k]>entmax
                                   entmax = ent[i,j,k]
                                   (index1, index2, index3) = (i,j,k)
                            end
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2], names(dataframe)[index3]
end

totalCorr = @time totalCorrelation(data)

print("total correlation: ", totalCorr)



# 8. returns PARTIAL INFORMATION DECOMPOSITION -------------DOES NOT TERMINATE!
function partialInformationDecomp(dataframe)
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
                            ent[i,j,k] = get_partial_information_decomposition(data_array[:,i],data_array[:,j],data_array[:,k])
                            if ent[i,j,k]>entmax
                                   entmax = ent[i,j,k]
                                   (index1, index2, index3) = (i,j,k)
                            end
                     end

              end
       end

       return entmax, names(dataframe)[index1], names(dataframe)[index2], names(dataframe)[index3]
end

partialInfoDecomp = @time partialInformationDecomp(data)

print("partial information decomposition: ", partialInfoDecomp)
