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

# Export our two super-useful functions
export find_primitive, refine_cell, standardize_cell,
        spg_get_dataset, get_spacegroup, get_symmetry,
        niggli_reduce, delaunay_reduce,
        ir_reciprocal_mesh

include("symmetrydb.jl")
include("standardize_cell.jl")
include("standardize_latt.jl")
include("ir_mesh.jl")


end #module LibSymspg
