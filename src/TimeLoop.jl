module TimeLoop

    using ArrayAllocators
    using NetCDF
    using Printf

    include("RunConfig.jl")
    include("IO.jl")
    include("Model.jl")
    include("GridBuffer.jl")
    include("Integrate.jl")

    export simulate

    function simulate(model, config)

        #creating the netCDF file to store our data
        IO.createNC(model, string(@__DIR__, "/../output/", config.filename), config)
    
        opPerc = config.steps/100 #number of steps correspondig to a percent of progress

        #Initialization. TODO replace by proper Initialization for n first steps required by model with lower order schemes
        buffer = GridBuffer.getTracerDiffInitBuffer(1000, model, config.gridSizes)

        for s = 1:config.steps

            #Showing the progress
            if s%opPerc == 0
                @printf "\033[5D%d%%" (s*100/config.steps) #\033[5D : escaped CSI control sequence to position the cursor
            end

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
    
        @printf "\n"

    end

end