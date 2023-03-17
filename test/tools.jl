using HDF5
import CSV 
using NBInclude

## convert notbook to script
function nb2jl(jl_name, nb_name=None)
    if nb_name == None
        nb_name = @__DIR__
    end
    nbexport(jl_name, nb_name)
end

## write csv 
function write_csv(save_path, df)
    CSV.write(save_path, df)
end

## write hdf5
function write_HDF5(file, key, value)
    try
        h5open(file, "r+") do f
            if haskey(f, key)
                delete_object(f, key)
            end
            f[key] = value
        end
        
    catch e

        h5open(file, "cw") do f
            f[key] = value
        end

    end
end

## write hdf5 with dict as inputs
function write_HDF5_dict(file, dicts)
    for (key, value) in dicts
        println(key)
        write_HDF5(file, key, value)
    end
end

function Retrieve_All_Ele(list)
    collect(Iterators.flatten(list))
end