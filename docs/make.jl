using SubgridTest
using Documenter

DocMeta.setdocmeta!(SubgridTest, :DocTestSetup, :(using SubgridTest); recursive=true)

makedocs(;
    modules=[SubgridTest],
    authors="fan.wang@tum.de",
    repo="https://github.com/FanWang0000/SubgridTest.jl/blob/{commit}{path}#{line}",
    sitename="SubgridTest.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
