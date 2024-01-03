module DeplacementFluide

using Printf

include("Model.jl")
include("IO.jl")
include("RunConfig.jl")
include("TimeLoop.jl")

export start

function start()
    model = Model.getLorenz()
    config = RunConfig.Config(10000, "a.nc", [3,3,3], :SSPRK3)

    @printf "Starting\n"

    @time TimeLoop.simulate(model, config)

end

end
