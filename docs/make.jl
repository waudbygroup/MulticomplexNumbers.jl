using Documenter
using MulticomplexNumbers

makedocs(
    sitename = "MulticomplexNumbers.jl",
    authors = "Chris Waudby",
    format = Documenter.HTML(),
    modules = [MulticomplexNumbers],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "User guide" => "userguide.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/waudbygroup/MulticomplexNumbers.jl.git"
)
