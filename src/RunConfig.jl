module RunConfig

export Config

struct Config
    steps::Integer
    filename::String
    gridSizes::Array{Integer}
end

end