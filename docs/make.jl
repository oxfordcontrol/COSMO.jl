using Documenter, DocumenterTools, COSMO

makedocs(
  sitename="COSMO.jl",
  format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://oxfordcontrol.github.io/COSMO.jl/stable/",),
  assets = ["assets/favicon.ico","assets/github_buttons.js"],
  authors = "Michael Garstka and contributors.",
  analytics = "UA-134239283-1",
  pages = [
        "Home" => "index.md",
        "User Guide" => "guide.md",
        "Method" => "method.md",
        "Examples" => "examples.md",
        "Citing COSMO" => "citing.md",
        "Contributing" => "contributing.md",
    ]
)

# deploydocs(
#     repo = "github.com/oxfordcontrol/COSMO.jl.git",
#     target = "build",
# )
