using SubgridTest
using Test

@testset "SubgridTest.jl" begin
    # Write your tests here.
    include("func.jl")
    include("func_LR.jl")
end
