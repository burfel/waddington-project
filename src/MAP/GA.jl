# -*- coding: utf-8 -*-
#=

The following program uses a genetic algorithm to approximate the mimimum action path (MAP).

Note: Call the functions by running f_g_functions_for_GA.jl

Args:
    -

Returns:
    -
=#

using Distributions, Plots, Bridge

# Defines a class Trajectory which contains
# - space::Array{Float64} -- the path through space
# - fitness::Float64 -- the fitess of the path depending on f and g
# - Nsize::Int64 -- the number of points that define the paths
type Trajectory
    space::Array{Float64}
    fitness::Float64
    Nsize::Int64
end


# Returns a straight line between two points given these two points and
# number of of points that define the path
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


# Returns a random path given a non_random_trajectory (see above)
# by adding a random value drawn from a normal distribution to each point
function random_traj(non_rand_traj)#,ID)
    #srand(ID)
    (Nsize,dim)=size(non_rand_traj)
    rand_traj=rand(Normal(0.0,0.5),Nsize,dim)
    #s = 1.
    #v = 0.
    for i=1:dim
        #B = sample(0:(1/(Nsize-1)):s, WienerBridge(s, v))
        #non_rand_traj[:,i]=non_rand_traj[:,i]+B.yy
        rand_traj[1,i]=0.0
        rand_traj[Nsize,i]=0.0
    end
    return non_rand_traj+rand_traj
end


# Initialises the population of trajectories for the GA
function pop_init(pop_size,pointA,pointB, Nsize,dim,f_func,g_func)
    # Initialises a new population array of data type Trajectory
    new_pop=Array{Trajectory}(pop_size)
    # Starting off with a straight line
    non_rand_traj=non_random_traj(pointA,pointB,Nsize,dim)
    # Computes the fitness of the straight line using f and g
    fit=fitness(non_rand_traj,f_func,g_func)
    # Creates first trajectory of initial population (a straight line)
    # TODO: vary the number of straight lines in initial poulation...
    new_pop[1]=Trajectory(non_rand_traj,fit,Nsize)
    # Creates the rest of the poulation by generating non random trajectories
    for i = 2:pop_size-1
        traj=random_traj(non_rand_traj)#,i)
        fit=fitness(traj,f_func,g_func)
        new_pop[i]=Trajectory(traj,fit,Nsize)
    end
    # smoothes the trajectories
    # TODO: try out different interpolation functions, eg cubic splines etc, see Interpolations.jl
    new_pop[pop_size]=smooth_traj(new_pop[1:pop_size-1],f_func,g_func)

    return new_pop
end


# Refines the path wrt f and g with a smoothing
# TODO: Specify smoothing...........
function refined_fitness(path,f_func,g_func,smoothing)
    (Nsize,dim)=size(path)
    new_path=zeros((Nsize-1)*smoothing+1,dim)

    for i=1:Nsize-1
        new_path[smoothing*(i-1)+1:smoothing*i,:]=non_random_traj(path[i,:],path[i+1,:],smoothing+1,dim)[1:smoothing,:]
    end
    new_path[(Nsize-1)*smoothing+1,:]=path[Nsize,:]

    return fitness(new_path,f_func,g_func)
end


# Computes the fitness, ie action of a path wrt to f and g
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


# Returns the child trajectory of two parent trajectories given f and g
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


# Returns the sorted population of a given size sorted_pop_size according to the indiviuals' fitness
# Note: sorted_pop_size must be <= pop_size
# TODO: Maybe add exception to test that requirement
# TODO: We could just use an implemented sort function...
function sorting(pop,sorted_pop_size)
    sorted_pop=Array{Trajectory}(sorted_pop_size)
    pop_size=size(pop)[1]
    l=zeros(sorted_pop_size)

    for i = 1:sorted_pop_size
        min=1000000000000000000000000000 #change this constant to the max possible
        k=0
        for j=1:pop_size
            if !(j in l) # if (j in l) or else noting
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


# Randomly selects two parents for the mating
function parents_selection(pop)
    p1=rand(1:size(pop)[1])
    p2=rand(1:size(pop)[1])

    parents=Array{Trajectory}(2)

    parents[1]=pop[p1]
    parents[2]=pop[p2]

    return parents
end


# Selects the new population from parents and children:
# pools the populations and selects the best ones
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


# Smoothes out the trajectories
# TODO: different methods of smooting, eg cubic splines etc
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


# Performs the GA
# TODO: to be discussed:
# - convergence
# - smoothing
# depending on initial population, mating (?)
# TODO: look into cubic splines and other interpolation methods to smooth the line integral
function GA(RepeatGA,Nb_gen,pop_size,pointA,pointB,Nsize,dim,f_func,g_func)
    repeatpop=pop_init(RepeatGA,pointA,pointB,Nsize,dim,f_func,g_func)
    r=0
    while r<RepeatGA
        pop=pop_init(pop_size,pointA,pointB,Nsize,dim,f_func,g_func)

        gen=0
        srand() # to set a random set to make results reproducible
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


    # Another GA with the best individuals of previous GA
    gen2=0
    while gen2<Nb_gen
        new_pop=Array{Trajectory}(size(repeatpop)[1])

        for i = 1:RepeatGA-1
            parents=parents_selection(repeatpop)
            new_pop[i]=mating(parents[1],parents[2],f_func,g_func)
        end
        new_pop[size(repeatpop)[1]]=smooth_traj(new_pop[1:size(repeatpop)[1]-1],f_func,g_func)
        repeatpop=selection(repeatpop,new_pop)
        gen2=gen2+1
    end

    return repeatpop
    #return smooth_traj(repeatpop,f_func,g_func)
end



#function mutation()

#end
# Another GA
# TODO: What are the differences? What are the results
function GA2(Nb_gen,pop_size,pointA,pointB,Nsize,dim,f_func,g_func)

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

    return pop
end



# Another GA
# What are the differences? What are the results?
function GA_mean(RepeatGA,Nb_gen,pop_size,pointA,pointB,Nsize,dim,f_func,g_func)
    #repeatpop=Array{Trajectory}(RepeatGA)
    repeatpop=pop_init(RepeatGA,pointA,pointB,Nsize,dim,f_func,g_func)
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
        r=r+1
        #repeatpop[r]=pop[1]
        repeatpop=selection(repeatpop,pop)
    end

    return smooth_traj(repeatpop,f_func,g_func)
end
