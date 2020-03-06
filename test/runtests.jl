using Test, Documenter, ModifiedHankelFunctionsOfOrderOneThird

println("Starting tests")

@time include("ModifiedHankelFunctionTests.jl")

println("Starting doctests")

doctest(ModifiedHankelFunctionsOfOrderOneThird)
