module Model

using Printf

export t_model
export t_variable
export get_lorenz

include("GridBuffer.jl")

struct t_variable #TODO add type to each variable ? (easier way to store vectors)
    name::String
    write::Bool #wether to write the variable to disk
    gridPos::Symbol #values taken on :face or :center
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
    variables[1] = t_variable("x(t)", true, :center)
    variables[2] = t_variable("y(t)", true, :center)
    variables[3] = t_variable("z(t)", true, :center)

    #update functions for the variables, in the same order

    RHS(vars) = permutedims([
        params[1] * (vars[2,:,:,:] .- vars[1,:,:,:]);;;; # σ[y(t)-x(t)]
        params[2] * vars[1,:,:,:] .- vars[2,:,:,:] .- vars[1,:,:,:] .* vars[3,:,:,:];;;; # ρx(t) - y(t) -x(t)z(t)
        vars[1,:,:,:] .* vars[2,:,:,:] .- params[3] * vars[3,:,:,:] # x(t)y(t) - βz(t)
    ], [4,1,2,3])

    m = t_model(params, 3, variables, 0.01, RHS)

    return m
end

function getTransport3DTraceurFlux()
    params = []

    variables = Array{t_variable}(undef, 4)
    variables[1] = t_variable("U", false, :face) #Vector type for variables ?
    variables[2] = t_variable("V", false, :face)
    variables[3] = t_variable("W", false, :face)
    variables[4] = t_variable("q", true, :center)

    #setting up update functions
    RHS(vars) = begin
        v1 = zeros(size(vars[1,:,:,:])...) #velocity field is constant, bad, should have constant fields as parameters I guess ?
        v2 = zeros(size(vars[1,:,:,:])...)
        v3 = zeros(size(vars[1,:,:,:])...)
        v4 = - GridBuffer.div(vars[4,:,:,:] .* vars[1,:,:,:], vars[4,:,:,:] .* vars[2,:,:,:], vars[4,:,:,:] .* vars[3,:,:,:]) # -∇·(Uq), U should be interpolated
        
        return permutedims(cat(v1,v2,v3,v4, dims=4), [4,1,2,3]) #ugly and unefficient TODO improve
    end

    m = t_model(params, 4, variables, 1, RHS)

    return m
end

function getTransport3DTraceurScalaire()
    params = []

    variables = Array{t_variable}(undef, 4)
    variables[1] = t_variable("U", false, :face) #Vector type for variables ?
    variables[2] = t_variable("V", false, :face)
    variables[3] = t_variable("W", false, :face)
    variables[4] = t_variable("q", true, :center)

    #setting up update functions
    RHS(vars) = begin
        v1 = zeros(size(vars[1,:,:,:])...) #velocity field is constant, bad, should have constant fields as parameters I guess ?
        v2 = zeros(size(vars[1,:,:,:])...)
        v3 = zeros(size(vars[1,:,:,:])...)
        
        gradq = GridBuffer.grad(vars[4,:,:,:])
        
        v4 = - (vars[1,:,:,:] .* gradq[1,:,:,:] .+ vars[2,:,:,:] .* gradq[2,:,:,:] .+ vars[3,:,:,:] .* gradq[3,:,:,:]) # -U·∇q, should be interpolated
        
        return permutedims(cat(v1,v2,v3,v4, dims=4), [4,1,2,3]) #ugly and unefficient TODO improve
    end

    m = t_model(params, 4, variables, 1, RHS)

    return m
end

end