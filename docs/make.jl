push!(LOAD_PATH,"../src/")

using Documenter, ModifiedHankelFunctionsOfOrderOneThird

DocMeta.setdocmeta!(ModifiedHankelFunctionsOfOrderOneThird, :DocTestSetup,
    :(using ModifiedHankelFunctionsOfOrderOneThird); recursive=true)

makedocs(
    modules=[ModifiedHankelFunctionsOfOrderOneThird],
    sitename="ModifiedHankelFunctionsOfOrderOneThird.jl",
    pages = [
        "Home" => "index.md",
        "API" => "pages/api.md",
        ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo="github.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl.git"
)
