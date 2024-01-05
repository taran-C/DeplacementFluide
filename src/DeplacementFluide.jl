module DeplacementFluide

using Printf

include("Model.jl")
include("IO.jl")
include("RunConfig.jl")
include("TimeLoop.jl")
include("GridBuffer.jl")

export start

function start()
    model = Model.getTransport3DTraceurFlux()
    config = RunConfig.Config(1000, "a.nc", [100,100,10], :SSPRK3, 100)

    @printf "Starting\n"

    @time TimeLoop.simulate(model, config, GridBuffer.getTracerDiffInitBuffer)

end

end
