module GridBuffer

using Printf

export Buffer
export interpolate

mutable struct Buffer
    size::Integer
    data::Array{Float64,5}
    currentIndex::Integer
    model
end

#Functions of an array q and an index i giving the interpolated value at q(i+1/2)
interpFuncs =   [(q, i)->q[i],
                (q, i)->(q[i]+q[i+1])/2.,
                (q, i)->(-q[i-1]+5*q[i]+2*q[i+1])/6.,
                (q, i)->(-q[i+1]+7*q[i]+7*q[i+1]-q[i+2])/12.,
                (q, i)->(2*q[i-2]-13*q[i-1]+47*q[i]+27*q[i+1]-3*q[i+2])/60.]

#computes q(i+1/2) with the maximum possible degree that we have a function for (TODO optimize)
#q : array along wich to interpolate
#i : index to get
#n : max degree that we want
function interpolate(q, i, n)
    deg = min(i, n, length(q) - i + 1, length(interpFuncs))
    return interpFuncs[deg](q, i)
end

#computes a FDA of dq at the point i
function fda(q, i)
    if i<length(q)
        return (q[i]+q[i+1])/2
    else 
        return q[i]
    end
end

function div3D(U, V, W) #to vectorize
    res = zeros(size(U)...)

    for i = 1:size(U)[1]-1, j = 1:size(V)[2]-1, k = 1:size(W)[3]-1
        res[i,j,k] = U[i,j,k] .- U[i+1,j,k] .+ V[i,j,k] .- V[i,j+1,k] .+ W[i,j,k] .- W[i,j,k+1]
    end
    #@printf "%s\n" res[2,2,2]
    return res
end

#Returns array of values initialized at one TODO function to initialize array from array of timesteps
function getOneBuffer(size::Integer, model, gridSizes::Array{Integer})
    data::Array{Float64,5} = ones(model.varNum,size,gridSizes...) #var t x y z
    return buffer = Buffer(size, data, 1, model)
end

#Returns array of values initialized at one TODO function to initialize array from array of timesteps
function getTracerDiffInitBuffer(size::Integer, model, gridSizes::Array{Integer})
    data = Array{Float64,5}(undef, model.varNum, size, gridSizes...)

    for v = 1:3, i = 1:gridSizes[1], j = 1:gridSizes[2], k = 1:gridSizes[3]
        data[v, 1, i,j,k] = sin(i)*sin(j)*sin(k)
    end

    #initial tracer distribution
    for i = 1:gridSizes[1], j = 1:gridSizes[2], k = 1:gridSizes[3]
        data[4, 1, i,j,k] = mod(floor(i/6),2) + mod(floor(j/6),2) + mod(floor(k/6),2)
    end

    return buffer = Buffer(size, data, 1, model)
end

end