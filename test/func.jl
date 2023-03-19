cd(@__DIR__)
# using SubgridTest
using DifferentialEquations, JLD2, Test

u0Two_LR = jldopen("./unittest_data.jld2", "r") do file
    read(file, "u0Two_LR")
end

exact_solTwoX_LR= jldopen("./unittest_data.jld2", "r") do file
    read(file, "solTwoX_LR")
end

exact_solTwoY_LR = jldopen("./unittest_data.jld2", "r") do file
    read(file, "solTwoY_LR")
end


u0One = jldopen("./unittest_data.jld2", "r") do file
    read(file, "u0One")
end


exact_solOneX = jldopen("./unittest_data.jld2", "r") do file
    read(file, "SolOneX")
end


struct VarSubgrid
    h::Real
    c::Real
end

# include("./tools.jl")

T0 = 0
Tmax = 10
K = 36 # so 10 degrees of longtitude per node
F = 10
J = 10

# ux0 = rand(K)
dt = 0.05

# uy0 = rand(K*J)
h = 1
c = 10
b = 10
# p = F, K, J, h, b, c
# args = p
psubgrid = VarSubgrid(h, c)
args = (K, F, J, b, psubgrid)
# args = (F, K, J, h, b, c)

# u0Two = [ux0; uy0]
argsOne = (F,K)
probOne = ODEProblem(Lorenz96One_shift!, u0One, (T0,Tmax), argsOne);
solOne = solve(probOne, dt=dt,saveat=dt);
solOneX = solOne[:,:]

probTwo = ODEProblem(Lorenz96Two_shift_LR!, u0Two_LR, (T0,Tmax), args);
solTwo = solve(probTwo, dt=dt, saveat=dt);

solTwo_t = solTwo.t
solTwoX = solTwo[1:K,:]
solTwoY = solTwo[K+1:end,:]
println(size(solOne))
println(size(exact_solOneX))
@test solOneX ≈ exact_solOneX
@test solTwoX≈exact_solTwoX_LR
@test solTwoY≈exact_solTwoY_LR