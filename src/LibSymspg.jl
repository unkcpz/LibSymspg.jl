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
export spg_find_primitive, spg_refine_cell, spg_get_dataset

mutable struct SpglibDB
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11, Cchar}
    hall_symbol::NTuple{17, Cchar}
    choice::NTuple{6, Cchar}
    transformation_matrix::NTuple{9, Float64}
    origin_shift::NTuple{3, Float64}
    n_operations::Cint
    rotations::Ptr{NTuple{9, Cint}}
    translations::Ptr{NTuple{3, Float64}}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{Tuple{7, Cchar}}
    equivalent_atoms::Ptr{Cint}
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{9, Float64}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3, Float64}}
    std_rotation_matrix::Ptr{NTuple{9, Float64}}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6, Cchar}
end

"""
spg_standardize_cell standardize cell

lattice are represented as row vector here
means each colume is a basis of lattice
positions are also represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function spg_standardize_cell(lattice::Array{Float64, 2},
                              positions::Array{Float64, 2},
                              types::Array{Int64, 1},
                              num_atom::Int64,
                              to_primitive::Int64,
                              no_idealize::Int64,
                              symprec::Float64=1e-5)
    # transpose to colomn vector used by spglib
    lattice = Array{Float64, 2}(lattice')
    positions = Array{Float64, 2}(positions')

    allocN = 4
    positions_ = zeros(Float64, 3, num_atom*allocN)
    types_ = zeros(Int64, num_atom*allocN)

    positions_[:, 1:num_atom] = positions[:, :]
    types_[1:num_atom] = types

    types_ = Base.cconvert(Array{Int32, 1}, types_)
    num_atom = Base.cconvert(Int32, num_atom)
    to_primitive = Base.cconvert(Int32, to_primitive)
    no_idealize = Base.cconvert(Int32, no_idealize)

    num_primitive_atom =
    ccall( (:spg_standardize_cell, libsymspg), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Int32, Int32, Float64 ),
           lattice, positions_, types_, num_atom, to_primitive, no_idealize, symprec )

    positions = positions_[:, 1:num_primitive_atom]
    types = types_[1:num_primitive_atom]

    # transpose back to row vectors
    lattice = Array{Float64, 2}(lattice')
    positions = Array{Float64, 2}(positions')

    return lattice, positions, Base.cconvert(Array{Int64, 1}, types), Base.cconvert(Int64,num_primitive_atom)
end

function spg_standardize_cell(lattice::Array{Float64, 2},
                              positions::Array{Float64, 2},
                              types::Array{Int64, 1},
                              num_atom::Int64,
                              to_primitive::Bool,
                              no_idealize::Bool,
                              symprec::Float64=1e-5)

    to_primitive ? to_primitive_ = 1 : to_primitive_ = 0
    no_idealize ? no_idealize_ = 1 : no_idealize_ = 0

    return spg_standardize_cell(lattice, positions, types, num_atom, to_primitive_, no_idealize_, symprec)
end

"""
spg_find_primitive reduce crystal to its primitive

lattice are represented as colume vector here
means each colume is a basis of lattice
while positions are represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function spg_find_primitive(lattice::Array{Float64, 2},
                            positions::Array{Float64, 2},
                            types::Array{Int64, 1},
                            num_atom::Int64,
                            symprec::Float64=1e-5)

    return spg_standardize_cell(lattice, positions, types, num_atom, 1, 0, symprec)
end

"""
spg_refine_cell refine crystal to its convetional

lattice are represented as colume vector here
means each colume is a basis of lattice
while positions are represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function spg_refine_cell(lattice::Array{Float64, 2},
                         positions::Array{Float64, 2},
                         types::Array{Int64, 1},
                         num_atom::Int64,
                         symprec::Float64=1e-5)

    return spg_standardize_cell(lattice, positions, types, num_atom, 0, 0, symprec)
end

function spg_get_dataset(lattice::Array{Float64, 2},
                         positions::Array{Float64, 2},
                         types::Array{Int64, 1},
                         num_atom::Int64,
                         symprec::Float64=1e-5)
    # transpose to colomn vector used by spglib
    lattice = Array{Float64, 2}(lattice')
    positions = Array{Float64, 2}(positions')

    types = Base.cconvert(Array{Int32, 1}, types)
    num_atom = Base.cconvert(Int32, num_atom)

    ptr_cdb =
    ccall( (:spg_get_dataset, libsymspg), Ptr{SpglibDB},
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    db = unsafe_load(ptr_cdb)
    println(db)
    return 1
end

end #module LibSymspg
