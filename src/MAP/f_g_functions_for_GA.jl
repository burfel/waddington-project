function f_test(vector)
    results=zeros(2)

    results[1]=vector[1]^4/(vector[1]^4+0.5^4)+0.5^4/(vector[2]^4+0.5^4)-vector[1]
    results[2]=vector[2]^4/(vector[2]^4+0.5^4)+0.5^4/(vector[1]^4+0.5^4)-vector[2]

    return results
end

function g_test(vector)
    results=zeros(2)

    results[1]=sqrt(abs(vector[1]^4/(vector[1]^4+0.5^4)+0.5^4/(vector[2]^4+0.5^4)-vector[1]))
    results[2]=sqrt(abs(vector[2]^4/(vector[2]^4+0.5^4)+0.5^4/(vector[1]^4+0.5^4)-vector[2]))

    return results
end

pointA=[0.0,3.0]
pointA=[0.5,1.5]
pointB=[1.0,1.0]
Nsize=40

pop=pop_init(10,pointA,pointB,Nsize,2,f_test,g_test)

function return_paths_from_pop(pop,id)
    results=zeros(size(pop[id].space)[1],size(pop[id].space)[2])

    for j=1:size(pop[id].space)[1]
        for k=1:size(pop[id].space)[2]
            results[j,k]=pop[id].space[j,k]
        end
    end
    return results
end

for i=1:size(pop)[1]
    paths=return_paths_from_pop(pop,i)
    plot(paths[:,1],paths[:,3])
end

pop=GA(20,100,100,pointA,pointB,Nsize,2,f_test,g_test)
paths=return_paths_from_pop(pop,1)
plot(paths[:,1],paths[:,2],xlims=[0.0,3.0],ylims=[0.0,3.0])




function f_3Dtest(vector)
    results=zeros(3)

    results[1]=vector[3]*vector[1]^4/(vector[1]^4+0.5^4)+0.5^4/(vector[2]^4+0.5^4)-vector[1]
    results[2]=vector[3]*vector[2]^4/(vector[2]^4+0.5^4)+0.5^4/(vector[1]^4+0.5^4)-vector[2]
    results[3]=-0.01*vector[3]

    return results
end

function g_3Dtest(vector)
    results=zeros(3)

    results[1]=sqrt(abs(vector[3]*vector[1]^4/(vector[1]^4+0.5^4)+0.5^4/(vector[2]^4+0.5^4)-vector[1]))
    results[2]=sqrt(abs(vector[3]*vector[2]^4/(vector[2]^4+0.5^4)+0.5^4/(vector[1]^4+0.5^4)-vector[2]))
    results[3]=sqrt(abs(-0.01*vector[3]))

    return results
end

pointA=[0.0,1.0,1.5]
pointB=[1.5,2.0,0.0]

pop=pop_init(10,pointA,pointB,Nsize,3,f_3Dtest,g_3Dtest)
pop=GA(20,100,100,pointA,pointB,Nsize,3,f_3Dtest,g_3Dtest)
paths=return_paths_from_pop(pop,9)
plotly()
plot(paths[:,1],paths[:,2],paths[:,3],xlims=[0.0,2.0],ylims=[0.0,2.0],zlims=[0.0,1.8])
