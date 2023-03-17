using Test, MyTestTwo # import packages that all tests need 

# maybe you load some data or define some constant that all tests need 

@testset "Test MyTestTwo.jl Basics" begin
    include("function_tests.jl") # include individual tests 
end 