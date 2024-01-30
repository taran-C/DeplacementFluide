module TimeLoop

    using ArrayAllocators
    using NetCDF
    using Printf
    using Dates

    include("RunConfig.jl")
    include("IO.jl")
    include("Model.jl")
    include("GridBuffer.jl")
    include("Integrate.jl")

    export simulate

    function simulate(model, config, bufferInit)

        #creating the netCDF file to store our data
        IO.createNC(model, string(@__DIR__, "/../output/", config.filename), config)

        #Initialization.
        buffer = bufferInit(config.buffSize, model, config.gridSizes)

        #setting up the progress display
        @printf "Starting...\n\n"
        start = Dates.Time(Dates.now())

        for s = 1:config.steps

            IO.showProgress(start, s, config.steps)

            Integrate.step!(config, model, buffer)


            if buffer.currentIndex == buffer.size -1
                # Write to NC file
                for v = 1:model.varNum
                    if model.variables[v].write
                        NetCDF.ncwrite(buffer.data[v,:,:,:,:], string(@__DIR__, "/../output/", config.filename), model.variables[v].name, start = [s - buffer.size + 2,1,1,1])
                    end
                end
                buffer.currentIndex = 0
            else
                buffer.currentIndex += 1
            end

        end

        #Writing what remains in the buffer
        for v = 1:model.varNum
            if model.variables[v].write
                NetCDF.ncwrite(buffer.data[v,:,:,:,:], string(@__DIR__, "/../output/", config.filename), model.variables[v].name, start = [config.steps-config.steps%buffer.size,1,1,1], count = [buffer.currentIndex,-1,-1,-1])
            end
        end
    
        @printf "Done !\n"

    end

end