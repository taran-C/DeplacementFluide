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

#Returns array of values initialized at one TODO function to initialize array from array of timesteps
function getOneBuffer(size::Integer, model, gridSizes::Array{Integer})
    data::Array{Float64,5} = ones(model.varNum,size,gridSizes...) #var t x y z
    return buffer = Buffer(size, data, 1, model)
end

#returns a certain value with interpolation if necessary TODO remove I guess, too specific
# function get(buffer, varNum::Integer, interpOrder::Integer, mode::Symbol, t::Integer, i::Integer, j::Integer, k::Integer)
#     if t<=0
#         t = t%buffer.size+buffer.size
#     end

#     if mode!=buffer.model.variables[varNum].gridPos
#         return interpolate(buffer.data[varNum, t, i, j, :], k, interpOrder)
#     else
#         return buffer.data[varNum, t, i, j, k]
#     end
# end

end