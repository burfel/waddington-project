
function Span_2D(boundaries, nb_bin, dim1, dim2)
    xspan=zeros(2,nb_bin)
    span_dim1=linspace(boundaries[1,dim1],boundaries[2,dim1],nb_bin)
    span_dim2=linspace(boundaries[1,dim2],boundaries[2,dim2],nb_bin)
    for bin=1:nb_bin
        xspan[1,bin]=span_dim1[bin]
        xspan[2,bin]=span_dim2[bin]
    end

    return xspan
end


function PFM_reduction_2D(end_state, boundaries, nb_bin, dim1, dim2)

    dim=size(end_state)[2]
    nb_sim=size(end_state)[1]

    Ps=zeros(Float64,nb_bin,nb_bin)
    total=0

    xspan=Span_2D(boundaries,nb_bin, dim1,dim2)

    for i = 1:nb_sim

        x1=1
        x2=1
        found1=false
        found2=false
        impossible=false

        while(!(found1 && found2))
            if (end_state[i,dim1]<xspan[1,1] || end_state[i,dim1]>=xspan[1,nb_bin] || end_state[i,dim2]<xspan[2,1] || end_state[i,dim2]>=xspan[2,nb_bin])
                impossible=true
                break;
            end


            if !(xspan[1,x1]<=end_state[i,dim1] && end_state[i,dim1]<xspan[1,x1+1]) && !found1
                x1=x1+1
            else
                found1=true
            end

            if !(xspan[2,x2]<=end_state[i,dim2] && end_state[i,dim2]<xspan[2,x2+1]) && !found2
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

    temp=(boundaries[2,dim1]-boundaries[1,dim1])*(boundaries[2,dim2]-boundaries[1,dim2])/(nb_bin*nb_bin)
    Ps=Ps/(total*temp)+0.00001

    return -log.(Ps)
end
