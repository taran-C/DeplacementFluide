module RunConfig

export Config

struct Config
    steps::Integer
    filename::String
    gridSizes::Array{Integer}
    integrator::Symbol
end

end