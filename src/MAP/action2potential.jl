## DEPRECATED -- TO BE DELETED!

# function that computes the action potential to its predecessor

#

# traj = Array{Trajectory}(1000)
traj = pop[1]
action = zeros(traj.Nsize, length(pop))
#plot(action)

for el in 1:length(pop)
    traj = pop[el]

    for elem in 1:traj.Nsize-1
        #action[elem] += traj[:fitness][elem]
        action[elem, el] = fitness(traj.space[elem:(elem+1),:], f, g)
    end
#    action
#    plot!(action)
end

#=
plot(actioProjection_2PCsn[:,1])
for i in length(pop)
    plot!(action[:,i])
end
=#

##-------------------------------------

traj = pop[1]
action = zeros(traj.Nsize)
#plot(action)

for elem in 1:traj.Nsize-1
        #action[elem] += traj[:fitness][elem]
        action[elem] = fitness(traj.space[elem:(elem+1),:], f, g)
end
#    action
#    plot!(action)



# plot the potential in 2D
plot(action)



# plot the potential in 3D (contour plot)
