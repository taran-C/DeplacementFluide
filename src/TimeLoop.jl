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

            showProgress(start, s, config.steps)

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

    function showProgress(start, s, steps)
        #Showing the progress
        updateEvery = 2
        if s%updateEvery == 0 || s==steps
            progress = s*100/steps
            progressInt = Int(floor(progress))

            elapsed = Dates.value(Dates.Millisecond(Dates.Time(Dates.now()) - start))
            minutesElapsed = elapsed÷60000
            secondsElapsed = (elapsed%60000)÷1000

            ETA = Int(floor(elapsed / progress * (100-progress)))
            minutesETA = ETA÷60000
            secondsETA = (ETA%60000)÷1000

            arrow = ">"
            if progress == 100
                arrow = ""
            end

            progressbar = "[" * "="^((progressInt)÷2) * arrow * " "^(50-progressInt÷2) * "]"
            @printf "\033[A\033[2K\033[G%2.2d%% %s Remaining %2.2d:%2.2d, elapsed %2.2d:%2.2d\n" progress progressbar minutesETA secondsETA minutesElapsed secondsElapsed #\033[G : escaped CSI control sequence to position the cursor \033[2K clears the line
        end
    end

end