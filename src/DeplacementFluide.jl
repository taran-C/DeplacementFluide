module DeplacementFluide

using Printf

include("Model.jl")
include("IO.jl")
include("RunConfig.jl")
include("TimeLoop.jl")

export start

function start()
    model = Model.getTransport3DTraceurFlux()
    config = RunConfig.Config(1000, "a.nc", [30,30,30], :SSPRK3)

    @printf "Starting\n"

    @time TimeLoop.simulate(model, config)

end

end
