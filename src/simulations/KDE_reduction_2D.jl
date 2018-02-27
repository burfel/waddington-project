

#we can discuss whether it is better to return the function KDE or the discretized matrix (in the format we use at the level of the code)
function KDE_reduction_2D(end_state, boundaries, nb_bin, dim1,dim2)

    dim=size(end_state)[2]
    nb_sim=size(end_state)[1]
    U_kde=zeros(Float64,nb_bin,nb_bin)

    #bandwidth
    hps=1.06*1000^(-1/5)

    #define the kernel density function from the end_state
    function KDE(x1,x2)
        sum = 0
        for i = 1:nb_sim
            #use the gaussian kernel (exp)
            sum=sum+exp(-((x1-end_state[i,dim1])^2+(x2-end_state[i,dim2])^2)/(2*(1.06*1000^(-1/5))^2))
        end

        #return the -log(Ps)
        sum=-log(sum/(1000*hps*sqrt(2*pi)))
        return sum
    end

    #define the space over which the function KDE is evaluated
    xspan=Span_2D(boundaries,nb_bin,dim1,dim2)

    #create the discretized matrix of KDE over the space previously defined
    for i = 1:nb_bin
        for j = 1:nb_bin
            U_kde[i,j]=KDE(xspan[1,i],xspan[2,j])
        end
    end

    #return the discretized matrix
    return U_kde

end
