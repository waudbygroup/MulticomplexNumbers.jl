using Documenter
using MulticomplexNumbers

makedocs(
    sitename = "MulticomplexNumbers.jl",
    authors = "Chris Waudby",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://waudbygroup.github.io/MulticomplexNumbers.jl/stable/",
    ),
    modules = [MulticomplexNumbers],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Background" => "background.md",
        "User Guide" => [
            "Creating Numbers" => "guide/creating.md",
            "Arithmetic Operations" => "guide/arithmetic.md",
            "Accessing Components" => "guide/components.md",
            "FFT Operations" => "guide/fft.md"
        ],
        "Applications" => [
            "NMR Spectroscopy" => "applications/nmr.md",
            "Numerical Differentiation" => "applications/differentiation.md"
        ],
        "Examples" => "examples.md",
        "API Reference" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/waudbygroup/MulticomplexNumbers.jl.git"
)
