push!(LOAD_PATH, string(@__DIR__, "/.."))
using GraphMatFun
using Documenter

DocMeta.setdocmeta!(GraphMatFun, :DocTestSetup, :(using GraphMatFun); recursive = true)

makedocs(;
    modules = [GraphMatFun],
    authors = "Elias Jarlebring <jarlebring@gmail.com> and contributors",
    repo = "https://github.com/matrixfunctions/GraphMatFun.jl/" *
           "blob/{commit}{path}#{line}",
    sitename = "GraphMatFun.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://matrixfunctions.github.io/GraphMatFun.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => ["graph.md", "generators.md", "optim.md", "code_gen.md"],
        "Index" => "docidx.md",
    ],
)

deploydocs(; repo = "github.com/matrixfunctions/GraphMatFun.jl", devbranch = "main")
