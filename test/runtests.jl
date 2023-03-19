using Pkg
# Pkg.add(path="https://github.com/FanWang0000/SubgridTest")
# Pkg.activate("..")

using SubgridTest
using Test

@testset "SubgridTest.jl" begin
    # Write your tests here.
    include("func.jl")
    include("func_LR.jl")
end
