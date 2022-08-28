using Documenter
using MulticomplexNumbers

makedocs(
    sitename = "MulticomplexNumbers",
    format = Documenter.HTML(),
    modules = [MulticomplexNumbers]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/waudbygroup/MulticomplexNumbers.jl.git"
)
