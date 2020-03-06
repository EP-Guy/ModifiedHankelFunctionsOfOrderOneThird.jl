using Test, Documenter, ModifiedHankelFunctionsOfOrderOneThird

println("Starting tests")

@time include("ModifiedHankelFunctionTests.jl")

println("Starting doctests")

DocMeta.setdocmeta!(ModifiedHankelFunctionsOfOrderOneThird, :DocTestSetup,
    :(using ModifiedHankelFunctionsOfOrderOneThird); recursive=true)

doctest(ModifiedHankelFunctionsOfOrderOneThird)
