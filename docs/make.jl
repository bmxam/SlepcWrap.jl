push!(LOAD_PATH, "../src/")

using SlepcWrap
using Documenter
using Literate

# Generate examples
example_src = joinpath(@__DIR__, "..", "example")
example_dir = joinpath(@__DIR__, "src", "example")
Sys.rm(example_dir; recursive=true, force=true)
Literate.markdown(joinpath(example_src, "helmholtz.jl"), example_dir; documenter=false, execute=false) # documenter = false to avoid Documenter to execute cells
Literate.markdown(joinpath(example_src, "helmholtz_fancy.jl"), example_dir; documenter=false, execute=false) # documenter = false to avoid Documenter to execute cells
Literate.markdown(joinpath(example_src, "complex.jl"), example_dir; documenter=false, execute=false) # documenter = false to avoid Documenter to execute cells
Literate.markdown(joinpath(example_src, "demo1.jl"), example_dir; documenter=false, execute=false) # documenter = false to avoid Documenter to execute cells

makedocs(;
    modules=[SlepcWrap],
    authors="bmxam",
    sitename="SlepcWrap.jl",
    clean=true, doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/bmxam/SlepcWrap.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => Any[
            "example/helmholtz.md",
            "example/demo1.md",
        ],
        "Fancy examples" => Any[
            "example/helmholtz_fancy.md",
            "example/complex.md",
        ],
        "API Reference" => Any[
            "api/init.md",
            "api/eps.md",
        ],
        "API fancy" => "api/fancy/fancy.md",
    ]
)

deploydocs(;
    repo="github.com/bmxam/SlepcWrap.jl.git",
    push_preview=true
)
