using Documenter, DocumenterTools, COSMO, Literate

@info "Building example problems..."

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
fix_suffix(filename) = replace(filename, ".jl" => ".md")
function postprocess(cont)
      """
      The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
      """ * cont
end
# find all example source files
exclude_files = ["chordal_decomposition.jl"; "chordal_decomposition_generate_data.jl"; "maxEigenvalue.jl"; "portfolio_optimization.jl"; "sum_of_squares.jl"];
example_path = joinpath(@__DIR__, "../examples/")
build_path =  joinpath(@__DIR__, "src", "examples/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)
filter!(x -> !in(x, exclude_files), files)

for file in files
      Literate.markdown(example_path * file, build_path; preprocess = fix_math_md, postprocess = postprocess, documenter = true, credit = true)
end

examples_nav = fix_suffix.("./examples/" .* files)
push!(examples_nav, "logistic_regression.md")

# find all other documentation source files that are build with Literate
example_path = joinpath(@__DIR__, "src", "literate/")
build_path =  joinpath(@__DIR__, "src", "literate", "build/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)
for file in files
      Literate.markdown(example_path * file, build_path; preprocess = fix_math_md, documenter = true, credit = true)
end


@info "Makeing documentation..."
makedocs(
  sitename="COSMO.jl",
  authors = "Michael Garstka and contributors.",
  format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://oxfordcontrol.github.io/COSMO.jl/stable/",
        assets = ["assets/favicon.ico","assets/github_buttons.js", "assets/cosmo.css"],
        analytics = "UA-134239283-1",
  ),
  pages = [
        "Home" => "index.md",
        "User Guide" => Any[
        "Getting Started" =>  "getting_started.md",
        "JuMP Interface" => "jump.md",
        "Linear System Solver" => "lin_solver.md",
        "Chordal Decomposition" => "decomposition.md",
        "Arbitrary Precision" => "./literate/build/arbitrary_precision.md"
        ],
        "Method" => "method.md",
        "Examples" => examples_nav,
        "Citing COSMO" => "citing.md",
        "Contributing" => "contributing.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/oxfordcontrol/COSMO.jl.git",
    target = "build")
