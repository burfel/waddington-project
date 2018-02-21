


function KDE_reduction_2D(end_state, boundaries, nb_bin, dim1,dim2)

    dim=size(end_state)[2]
    nb_sim=size(end_state)[1]
    U_kde=zeros(Float64,nb_bin,nb_bin)

    hps=1.06*1000^(-1/5)

    function KDE(x1,x2)
        sum = 0
        for i = 1:nb_sim
            sum=sum+exp(-((x1-end_state[i,dim1])^2+(x2-end_state[i,dim2])^2)/(2*(1.06*1000^(-1/5))^2))
        end
        sum=-log(sum/(1000*hps*sqrt(2*pi)))
        return sum
    end

    xspan=Span_2D(boundaries,nb_bin,dim1,dim2)

    for i = 1:nb_bin
        for j = 1:nb_bin
            U_kde[i,j]=KDE(xspan[1,i],xspan[2,j])
        end
    end

    return U_kde

end
