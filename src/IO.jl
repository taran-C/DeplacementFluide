module IO

    using NetCDF
    using Printf

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
end