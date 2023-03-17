module SubgridTest
    # using DifferentialEquations
    include("L96func.jl") 
    include("tools.jl") 
    include("tools_stats.jl") 
# Write your package code here.
    export KS_test,
           DiscreteDensity, 
           poly_fit,
           nb2jl, 
           write_csv, 
           write_HDF5, 
           write_HDF5_dict, 
           Retrieve_All_Ele,
           lorenz96!, 
           Lorenz96One_shift!, 
           Lorenz96Two_shift_LR!, 
           Lorenz96Two_shift!, 
           Lorenz96Three_shift!,
           Lorenz96Two_polyB!, 
           Lorenz96Two_polyBk!
    function __init__() # OPTIONAL: this special function is always executed when the module is loaded 
        nothing 
    end
end
