module RunConfig

export Config

struct Config
    steps::Integer
    filename::String
    gridSizes::Array{Integer}
    integrator::Symbol
    buffSize::Integer
    interpDeg::Integer
end

end