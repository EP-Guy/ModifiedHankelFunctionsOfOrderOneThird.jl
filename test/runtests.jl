using Test, Documenter, ModifiedHankelFunctionsOfOrderOneThird

println("\nStarting tests")

@time include("ModifiedHankelFunctionTests.jl")

println("\nStarting doctests")

DocMeta.setdocmeta!(ModifiedHankelFunctionsOfOrderOneThird, :DocTestSetup,
    :(using ModifiedHankelFunctionsOfOrderOneThird); recursive=true)

doctest(ModifiedHankelFunctionsOfOrderOneThird)
