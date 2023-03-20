cd(@__DIR__)
# using SubgridTest
using DifferentialEquations, JLD2, Test

# jldsave("unittest_data_LR.jld2"; BoneX=BoneX, BkX=BkX, u0_BLR=u0, density_BkX=density_BkX, density_BoneX)


exact_BoneX = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "BoneX")
end

exact_BkX= jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "BkX")
end

exact_density_BkX = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "density_BkX")
end


u0_BLR = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "u0_BLR")
end


exact_density_BoneX = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "density_BoneX")
end

B_poly_K = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "B_poly_K")
end

B_poly_one = jldopen("./unittest_data_LR.jld2", "r") do file
    read(file, "B_poly_one")
end

struct VarSubgrid
    h::Real
    c::Real
end

T0 = 0 
# end time
# Tmax = 100
Tmax = 10
K = 36 # so 10 degrees of longtitude per node
F = 10
# ux0 = rand(K)
dt = 0.05
sdt=dt
J = 10
# uy0 = rand(K*J)
h = 1
c = 10
b = 10
psubgrid = VarSubgrid(h, c)

argsBone = (K, F, J, b, psubgrid, B_poly_one)



argsBone = (K, F, J, b, psubgrid, B_poly_one)
argsBk = (K, F, J, b, psubgrid, B_poly_K)
# u0 = [ux0; uy0];
probTwoBone = ODEProblem(Lorenz96Two_polyB!, u0_BLR, (T0,Tmax), argsBone);
solTwoBone = solve(probTwoBone, dt=dt,saveat=dt);

# println(size(u0_BLR))
probTwoBk = ODEProblem(Lorenz96Two_polyBk!, u0_BLR, (T0,Tmax), argsBk);
solTwoBk = solve(probTwoBk, dt=dt, saveat=dt);

BoneX = solTwoBone[:, :]
BkX = solTwoBk[:,:]
# TwoB_estX = solTwoB_est[1:K,:];
println(exact_BoneX[end,end])
println(BoneX[end,end])
@testset "solution of Bx" begin
@test BoneX ≈ exact_BoneX
end 

@testset "solution of BOne" begin
@test BkX≈exact_BkX
end


bmax = maximum([BoneX[:]; BkX[:]])
bmin = minimum([BoneX[:]; BkX[:]])
bins = LinRange(bmin, bmax, 100);
# density_TwoB_estX = DiscreteDensity(TwoB_estX[:], bins);
density_BoneX = DiscreteDensity(BoneX[:], bins);
density_BkX = DiscreteDensity( BkX[:], bins);

@testset "empirical density estimation of BOne" begin
    @test density_BoneX ≈ exact_density_BoneX
    end

@testset "empirical density estimation of Bk" begin
@test density_BkX ≈ exact_density_BkX
end