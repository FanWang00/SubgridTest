module SubgridTest
    # using Pkg
    # Pkg.add("Pkg")
# Write your package code here.
    # using DifferentialEquations
    include("L96func.jl") 
    # include("tools.jl") 
    include("tools_stats.jl") 
# Write your package code here.
    export KS_test,
           DiscreteDensity, 
           poly_fit,

           Lorenz96One_shift!, 
           Lorenz96Two_shift_LR!, 
           Lorenz96Two_shift!, 
           
           Lorenz96Two_polyB!, 
           Lorenz96Two_polyBk!
    function __init__() # OPTIONAL: this special function is always executed when the module is loaded 
        nothing 
    end
end