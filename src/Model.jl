module Model

export t_model
export t_variable
export get_lorenz

include("GridBuffer.jl")

struct t_variable
    name::String
    write::Bool #wether to write the variable to disk
    gridPos::Symbol #values taken on :face or :center
end

struct t_model
    params::Array{AbstractFloat}
    varNum::Integer
    variables::Array{t_variable}
    deltaT::AbstractFloat
    updateFuncs::Array{Function}
end

function getLorenz()
    params = Array{Float64}(undef, 3)
    params[1] = 10.
    params[2] = 28.
    params[3] = 8. / 3.

    variables = Array{t_variable}(undef, 3)
    variables[1] = t_variable("x(t)", true, :center)
    variables[2] = t_variable("y(t)", true, :center)
    variables[3] = t_variable("z(t)", true, :center)

    #update functions for the variables, in the same order
    functions = Array{Function}(undef, 3)
    f1(model::t_model, buffer, i::Integer, j::Integer, k::Integer) = model.deltaT * model.params[1] * (GridBuffer.get(buffer, 2, 1, :center, buffer.currentIndex, i,j,k) - GridBuffer.get(buffer, 1, 1, :center, buffer.currentIndex, i,j,k)) + GridBuffer.get(buffer, 1, 1, :center, buffer.currentIndex, i,j,k)
    f2(model::t_model, buffer, i::Integer, j::Integer, k::Integer) = model.params[2] * model.deltaT * GridBuffer.get(buffer, 1, 1, :center, buffer.currentIndex, i,j,k) + GridBuffer.get(buffer, 2, 1, :center, buffer.currentIndex, i,j,k) * (1. - model.deltaT) - model.deltaT * GridBuffer.get(buffer, 1, 1, :center, buffer.currentIndex, i,j,k) * GridBuffer.get(buffer, 3, 1, :center, buffer.currentIndex, i,j,k)
    f3(model::t_model, buffer, i::Integer, j::Integer, k::Integer) = model.deltaT * GridBuffer.get(buffer, 1, 1, :center, buffer.currentIndex, i,j,k) * GridBuffer.get(buffer, 2, 1, :center, buffer.currentIndex, i,j,k) + GridBuffer.get(buffer, 3, 1, :center, buffer.currentIndex, i,j,k) * (1. - model.deltaT * model.params[3])
    functions[1] = f1
    functions[2] = f2
    functions[3] = f3

    m = t_model(params, 3, variables, 0.01, functions)

    return m
end

function getTransport3DTraceurFlux()
    params = []

    variables = Array{t_variable}(undef, 4)
    variables[1] = t_variable("U", true, :face)
    variables[2] = t_variable("V", true, :face)
    variables[3] = t_variable("W", true, :face)
    variables[4] = t_variable("q", :diagnostic, :center)

    #setting up update functions
    id(varNum) = (model::t_model, buffer, i::Integer, j::Integer, k::Integer) -> buffer[varNum, buffer.currentIndex, i,j,k]

    fq(model::t_model, buffer, i::Integer, j::Integer, k::Integer) = begin
        #interpolating the speed field to the center of the mesh volumes
        Uc = GridBuffer.interpolate(buffer.data[1, buffer.currentIndex, :, j, k], i, 3)
        Vc = GridBuffer.interpolate(buffer.data[2, buffer.currentIndex, i, :, k], j, 3)
        Wc = GridBuffer.interpolate(buffer.data[3, buffer.currentIndex, i, j, :], k, 3)

        #Fq = - (GridBuffer.fda(buffer.data[4, buffer.currentIndex]))

        q1 = GridBuffer.get(buffer, 4, 1, :center, buffer.currentIndex, i,j,k) + model.deltaT * 2
    end

    functions = Array{Function}(undef, 4)
    functions[1] = id(1)
    functions[2] = id(2)
    functions[3] = id(3)
    functions[4] = fq

end

end