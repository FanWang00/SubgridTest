using SubgridTest00
using Documenter

DocMeta.setdocmeta!(SubgridTest00, :DocTestSetup, :(using SubgridTest00); recursive=true)

makedocs(;
    modules=[SubgridTest00],
    authors="fan.wang@tum.de",
    repo="https://github.com/FanWang0000/SubgridTest00.jl/blob/{commit}{path}#{line}",
    sitename="SubgridTest00.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
