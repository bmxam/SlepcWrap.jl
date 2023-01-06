using SlepcWrap
using Test


# from :
# https://discourse.julialang.org/t/what-general-purpose-commands-do-you-usually-end-up-adding-to-your-projects/4889
@generated function compare_struct(x, y)
    if !isempty(fieldnames(x)) && x == y
        mapreduce(n -> :(x.$n == y.$n), (a, b) -> :($a && $b), fieldnames(x))
    else
        :(x == y)
    end
end


@testset "SlepcWrap.jl" begin
    include("./helmholtz.jl")
    include("./helmholtz_fancy.jl")
    include("./demo1.jl")
end