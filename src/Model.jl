module Model

export t_model
export t_variable
export get_lorenz

include("GridBuffer.jl")

struct t_variable
    name::String
    write::Bool #wether to write the variable to disk
    gridPos::Symbol #values taken on :face or :center
    type::Symbol
end

struct t_model
    params::Array{AbstractFloat}
    varNum::Integer
    variables::Array{t_variable}
    deltaT::AbstractFloat
    RHS::Function
end

function getLorenz()
    params = Array{Float64}(undef, 3)
    params[1] = 10. #sigma
    params[2] = 28. #rho
    params[3] = 8. / 3. #beta

    variables = Array{t_variable}(undef, 3)
    variables[1] = t_variable("x(t)", true, :center, :prognostic)
    variables[2] = t_variable("y(t)", true, :center, :prognostic)
    variables[3] = t_variable("z(t)", true, :center, :prognostic)

    #update functions for the variables, in the same order
    functions = Array{Function}(undef, 3)

    RHS(vars) = [params[1] * (vars[2] .- vars[1]), # σ[y(t)-x(t)]
                params[2] * vars[1] .- vars[2] .- vars[1] .* vars[3], # ρx(t) - y(t) -x(t)z(t)
                vars[1] .* vars[2] .- params[3] * vars[3]] # x(t)y(t) - βz(t)

    m = t_model(params, 3, variables, 0.01, RHS)

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