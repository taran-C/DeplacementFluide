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

#computes a FDA of dq along dimension dim, TODO upstream
function dq(q, dim)
    s = size(q)
    res = zeros(s...)

    if dim == 1
        res[2:s[1]-1,:,:] = (q[3:s[1],:,:]-q[1:s[1]-2,:,:])/2.
    elseif dim == 2
        res[:,2:s[2]-1,:] = (q[:,3:s[2],:]-q[:,1:s[2]-2,:])/2.
    elseif dim == 3
        res[:,:,2:s[3]-1] = (q[:,:,3:s[3]]-q[:,:,1:s[3]-2])/2.
    end

    # FDA for the second-order derivative
    # if dim == 1
    #     res[2:s[1]-1,:,:] = (q[3:s[1],:,:] .+ q[1:s[1]-2,:,:] .- 2. *q[2:s[1]-1,:,:])
    # elseif dim == 2
    #     res[:,2:s[2]-1,:] = (q[:,3:s[2],:] .+ q[:,1:s[2]-2,:] .- 2. *q[:,2:s[2]-1,:])
    # elseif dim == 3
    #     res[:,:,2:s[3]-1] = (q[:,:,3:s[3]] .+ q[:,:,1:s[3]-2] .- 2. *q[:,:,2:s[3]-1])
    # end

    return res
end

function div(U, V, W)
    return dq(U, 1) + dq(V, 2) + dq(W, 3)
end

function grad(q)
    return permutedims([dq(q,1);;;;dq(q,2);;;;dq(q,3)], [4,1,2,3])
end

function cross(u1,v1,w1,u2,v2,w2)
    return [v1.*w2 .- w1.*v2;;;; w1.*u2 .- u1.*w2;;;; u1.*v2 .- v1.*u2]
end

#Returns array of values initialized at one TODO function to initialize array from array of timesteps
function getOneBuffer(size::Integer, model, gridSizes::Array{Integer})
    data::Array{Float64,5} = ones(model.varNum,size,gridSizes...) #var t x y z
    return buffer = Buffer(size, data, 1, model)
end

#Returns array of values initialized at one TODO function to initialize array from array of timesteps
function getTracerDiffInitBuffer(size::Integer, model, gridSizes::Array{Integer})
    data = Array{Float64,5}(undef, model.varNum, size, gridSizes...)

    #velocity field
    psiA = Array{Float64,3}(undef, gridSizes...)
    psiB = Array{Float64,3}(undef, gridSizes...)

    for i = 1:gridSizes[1], j = 1:gridSizes[2], k = 1:gridSizes[3]
        psiA[i,j,k] = sin(pi*i/gridSizes[1])*sin(pi*j/gridSizes[2])
        psiB[i,j,k] = sin(pi*i/gridSizes[1])*sin(pi*k/gridSizes[3])
    end

    gradA = grad(psiA)
    gradB = grad(psiB)

    one = ones(gridSizes...)
    zero = zeros(gridSizes...)

    Ua = cross(zero,zero,one, gradA[1,:,:,:], gradA[2,:,:,:], gradA[3,:,:,:])
    Ub = cross(zero,one,zero, gradB[1,:,:,:], gradB[2,:,:,:], gradB[3,:,:,:])

    v = Ua + Ub

    v = permutedims(v, [4,1,2,3])

    data[1,1,:,:,:] = v[1,:,:,:]
    data[2,1,:,:,:] = v[2,:,:,:]
    data[3,1,:,:,:] = v[3,:,:,:]

    #initial tracer distribution
    for i = 1:gridSizes[1], j = 1:gridSizes[2], k = 1:gridSizes[3]
        data[4, 1, i,j,k] = mod(floor(i/6),2) + mod(floor(j/6),2) + mod(floor(k/6),2)
    end

    return buffer = Buffer(size, data, 1, model)
end

end