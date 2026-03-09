using Documenter
using SparseSensors

makedocs(
    sitename = "SparseSensors.jl",
    modules = [SparseSensors],
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/samtalki/SparseSensors.jl.git",
    devbranch = "main",
)
