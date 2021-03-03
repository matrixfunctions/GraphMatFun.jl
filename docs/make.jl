using GraphMatFun
using Documenter

DocMeta.setdocmeta!(GraphMatFun, :DocTestSetup, :(using GraphMatFun); recursive=true)

makedocs(;
    modules=[GraphMatFun],
    authors="Elias Jarlebring <jarlebring@gmail.com> and contributors",
    repo="https://github.com/jarlebring/GraphMatFun.jl/blob/{commit}{path}#{line}",
    sitename="GraphMatFun.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jarlebring.github.io/GraphMatFun.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jarlebring/GraphMatFun.jl",
)
