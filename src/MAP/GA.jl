using Distributions
using Plots


type Trajectory
    space::Array{Float64}
    fitness::Float64
    Nsize::Int64
end



function non_random_traj(pointA, pointB,Nsize,dim)
    non_rand_traj=zeros(Nsize,dim)
    for i = 1:dim
        xi=linspace(pointA[i],pointB[i],Nsize)
        for j = 1:Nsize
            non_rand_traj[j,i]=xi[j]
        end
    end
    return non_rand_traj
end



function random_traj(non_rand_traj)#,ID)
    #srand(ID)
    (Nsize,dim)=size(non_rand_traj)
    rand_traj=rand(Normal(0.0,0.5),Nsize,dim)
    for i=1:dim
        rand_traj[1,i]=0.0
        rand_traj[Nsize,i]=0.0
    end
    return non_rand_traj+rand_traj
end



function pop_init(pop_size,pointA,pointB, Nsize,dim,f_func,g_func)
    new_pop=Array{Trajectory}(pop_size)
    non_rand_traj=non_random_traj(pointA,pointB,Nsize,dim)

    fit=fitness(non_rand_traj,f_func,g_func)
    new_pop[1]=Trajectory(non_rand_traj,fit,Nsize)

    for i = 2:pop_size-1
        traj=random_traj(non_rand_traj)#,i)
        fit=fitness(traj,f_func,g_func)
        new_pop[i]=Trajectory(traj,fit,Nsize)
    end
    new_pop[pop_size]=smooth_traj(new_pop[1:pop_size-1],f_func,g_func)

    return new_pop
end


#fitness=action
function fitness(path,f_func,g_func)
    fitness = 0.0
    (Nsize,dim)=size(path)
    ds=zeros(Nsize-1,dim)
    f=zeros(Nsize-1,dim)
    g_adapted=zeros(Nsize-1,dim)
    for i=1:Nsize-1
        for j=1:dim
            ds[i,j]=path[i+1,j]-path[i,j]

            if (g_func(path[i,:])[j]==0)
                g_adapted[i,j]=0
            else
                g_adapted[i,j]=1/g_func(path[i,:])[j]^2
            end
        end
        f[i,:]=f_func(path[i,:])
    end
    for i=1:Nsize-1
        temp1=0.0
        temp2=0.0
        temp3=0.0
        ds_length=0.0
        for j=1:dim
            temp1=temp1+g_adapted[i,j]*f[i,j]^2
            temp2=temp2+g_adapted[i,j]*ds[i,j]^2
            temp3=temp3+g_adapted[i,j]*ds[i,j]*f[i,j]
            ds_length=ds_length+ds[i,j]^2
        end
        temp1=sqrt(temp1)
        temp1=sqrt(temp2)
        ds_length=sqrt(ds_length)
        if temp2==0
            fitness=fitness+temp1
        elseif !(temp1==0)
            fitness=fitness+abs(temp1*(1-temp3/(temp1*temp2))*ds_length)
        end
    end
    return 2*fitness
end





function mating(parent1,parent2,f_func,g_func)
    (Nsize,dim)=size(parent1.space)
    new_traj=zeros(Nsize,dim)
    cut=rand(1:Nsize)

    for i = 1:Nsize
        for j=1:dim
            if i<cut
                new_traj[i,j]=parent1.space[i,j]
            else
                new_traj[i,j]=parent2.space[i,j]
            end
        end
    end

    child = Trajectory(new_traj,fitness(new_traj,f_func,g_func),Nsize)
    return child
end


#new_pop_size <= pop_size
function sorting(pop,sorted_pop_size)
    sorted_pop=Array{Trajectory}(sorted_pop_size)
    pop_size=size(pop)[1]
    l=zeros(sorted_pop_size)

    for i = 1:sorted_pop_size
        min=1000000000000000000000000000 #change this constant to the max possible
        k=0
        for j=1:pop_size
            if !(j in l)
                if pop[j].fitness<min
                    min=pop[j].fitness
                    k=j
                end
            end
        end
        l[i]=k
        #if min==0
        #    for u = 1:sorted_pop_size
        #        sorted_pop[u]=pop[k]
        #    end
        #    break
        #else
        sorted_pop[i]=pop[k]
        #end
    end
    return sorted_pop
end



function parents_selection(pop)
    p1=rand(1:size(pop)[1])
    p2=rand(1:size(pop)[1])

    parents=Array{Trajectory}(2)

    parents[1]=pop[p1]
    parents[2]=pop[p2]

    return parents
end



function selection(parents,children)
    parents_size=size(parents)[1]
    children_size=size(children)[1]
    total=Array{Trajectory}(parents_size+children_size)
    for i = 1:(parents_size+children_size)
        if i <=parents_size
            total[i]=parents[i]
        else
            total[i]=children[i-parents_size]
        end
    end
    total=sorting(total,parents_size)
    return total
end




function smooth_traj(pop,f_func,g_func)
    (Nsize,dim)=size(pop[1].space)
    size_pop=size(pop)[1]
    mean_traj=zeros(Nsize,dim)

    for i=1:Nsize
        for j=1:dim
            for u=1:size_pop
                mean_traj[i,j]=mean_traj[i,j]+pop[u].space[i,j]
            end
            mean_traj[i,j]=mean_traj[i,j]/size_pop
        end
    end

    new_traj=Trajectory(mean_traj,fitness(mean_traj,f_func,g_func),Nsize)
    return new_traj
end

function GA(RepeatGA,Nb_gen,pop_size,pointA,pointB,Nsize,dim,f_func,g_func)
    repeatpop=pop_init(pop_size,pointA,pointB,Nsize,dim,f_func,g_func)
    r=0
    while r<RepeatGA
        pop=pop_init(pop_size,pointA,pointB,Nsize,dim,f_func,g_func)

        gen=0
        srand()
        while gen<Nb_gen
            new_pop=Array{Trajectory}(pop_size)

            for i = 1:pop_size-1
                parents=parents_selection(pop)
                new_pop[i]=mating(parents[1],parents[2],f_func,g_func)
            end
            new_pop[pop_size]=smooth_traj(new_pop[1:pop_size-1],f_func,g_func)

            pop=selection(pop,new_pop)
            gen=gen+1
        end
        repeatpop=selection(repeatpop,pop)
        r=r+1
    end

    gen2=0
    while gen2<Nb_gen
        new_pop=Array{Trajectory}(size(repeatpop)[1])

        for i = 1:pop_size-1
            parents=parents_selection(repeatpop)
            new_pop[i]=mating(parents[1],parents[2],f_func,g_func)
        end
        new_pop[size(repeatpop)[1]]=smooth_traj(new_pop[1:size(repeatpop)[1]-1],f_func,g_func)
        repeatpop=selection(repeatpop,new_pop)
        gen2=gen2+1
    end

    return repeatpop
end
#function mutation()

#end
