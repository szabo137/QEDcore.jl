using QEDcore
using Documenter

DocMeta.setdocmeta!(QEDcore, :DocTestSetup, :(using QEDcore); recursive=true)

makedocs(;
    modules=[QEDcore],
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de>",
    sitename="QEDcore.jl",
    format=Documenter.HTML(;
        canonical="https://QEDjl-project.github.io/QEDcore.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/QEDjl-project/QEDcore.jl", devbranch="main")
