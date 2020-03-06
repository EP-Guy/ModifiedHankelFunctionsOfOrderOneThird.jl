push!(LOAD_PATH,"../src/")

using Documenter, ModifiedHankelFunctionsOfOrderOneThird

makedocs(
    modules=[ModifiedHankelFunctionsOfOrderOneThird],
    sitename="ModifiedHankelFunctionsOfOrderOneThird.jl",
    pages = [
        "Home" => "index.md",
        "API" => "pages/api.md",
        ],
    )

deploydocs(
    repo="github.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl.git"
)
