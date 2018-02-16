
function Span(boundaries, nb_bin, dim)
    xspan=zeros(dim,nb_bin)
    for j = 1:dim
        span=linspace(boundaries[1,j],boundaries[2,j],nb_bin)
        for bin=1:nb_bin
            xspan[j,bin]=span[bin]
        end
    end

    return xspan
end


function Prob_Flux_Method(end_state, boundaries, tuple_bindim)

    dim=size(end_state)[2]
    nb_sim=size(end_state)[1]

    nb_bin=tuple_bindim[1]

    Ps=zeros(Float64,tuple_bindim)
    total=0

    xspan=Span(boundaries,nb_bin,dim)

    for i = 1:nb_sim

        x=ones(Int64,dim)
        found=falses(dim)
        impossible=false

        for j=1:dim
            if (end_state[i,j]<xspan[j,1] || end_state[i,j]>=xspan[j,nb_bin])
                impossible=true
                break;
            end
        end

        if !impossible
            while(!all(found))
                for j=1:dim
                    if !(xspan[j,x[j]]<=end_state[i,j] && end_state[i,j]<xspan[j,x[j]+1]) && !found[j]
                        x[j]=x[j]+1
                    else
                        found[j]=true
                    end
                end
            end

            if dim==1
                Ps[x[1]]=Ps[x[1]]+1
            elseif dim==2
                Ps[x[1],x[2]]=Ps[x[1],x[2]]+1
            elseif dim==3
                Ps[x[1],x[2],x[3]]=Ps[x[1],x[2],x[3]]+1
            elseif dim==4
                Ps[x[1],x[2],x[3],x[4]]=Ps[x[1],x[2],x[3],x[4]]+1
            elseif dim==5
                Ps[x[1],x[2],x[3],x[4],x[5]]=Ps[x[1],x[2],x[3],x[4],x[5]]+1
            end

            total=total+1
        end
    end

    temp=1
    for j=1:dim
        temp=temp*(boundaries[2,j]-boundaries[1,j])/nb_bin
    end

    Ps=Ps/(total*temp)+0.00001

    return -log.(Ps)
end
