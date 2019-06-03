module LibSymspg
using Libdl

# Load in `deps.jl`, complaining if it does not exist
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    error("LibSymspg not installed properly, run Pkg.build(\"LibSymspg\"), restart Julia and try again")
end
include(depsjl_path)

# Module initialization function
function __init__()
    # Always check your dependencies from `deps.jl`
    check_deps()
end

export char2Str, rotsFromTuple, transFromTuple

function char2Str(charTuple::Tuple{Vararg{UInt8}})
    cs = collect(charTuple)
    idx = 1
    for i in eachindex(cs)
        if cs[i] == 0
            break
        end
        idx = i
    end
    cs = Char.(cs[1:idx])

    return String(cs)
end

function rotsFromTuple(rotsTuple::Array{NTuple{9,Int32},1}, nop::Integer)
    r = Array{Int64,3}(undef, 3, 3, nop)
    for i in 1:nop
        r[:,:,i] = reshape([Base.convert(Int64, e) for e in rotsTuple[i]], 3, 3)
    end
    return r
end

function transFromTuple(transTuple::Array{NTuple{3,Float64}}, nop::Integer)
    t = Array{Float64,2}(undef, 3, nop)
    for i in 1:nop
        t[:,i] = [e for e in transTuple[i]]
    end
    return t
end

include("version.jl")
include("symmetry-api.jl")
include("spacegroup-api.jl")
include("cell-reduce-api.jl")
include("latt-reduce-api.jl")
include("ir-mesh-api.jl")


end #module LibSymspg
