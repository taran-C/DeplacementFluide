module IO

    using NetCDF
    using Printf
    using Dates

    export createNC

    function createNC(model, filename, config)
        if isfile(filename)
            rm(filename)
        end

        for i = 1:model.varNum
            if model.variables[i].write
                NetCDF.nccreate(filename, model.variables[i].name, "t", -1, "x", config.gridSizes[1], "y", config.gridSizes[2], "z", config.gridSizes[3])
            end
        end
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