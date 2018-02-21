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

pointB=[1.0,1.0]
pointA=[0.5,1.5]
Nsize=20

pop=pop_init(10,pointA,pointB,Nsize,2,f_test,g_test)
