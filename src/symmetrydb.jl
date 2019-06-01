export get_spacegroup

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

function spg_get_dataset(lattice::Array{Float64, 2},
                         positions::Array{Float64, 2},
                         types::Array{Int64, 1},
                         num_atom::Int64,
                         symprec::Float64=1e-5)

    @assert size(positions)[1] == size(types)[1] == num_atom

    # transpose to colomn vector used by spglib
    positions = Array{Float64, 2}(positions')

    types = Base.cconvert(Array{Int32, 1}, types)
    num_atom = Base.cconvert(Int32, num_atom)

    ptr_cdb =
    ccall( (:spg_get_dataset, libsymspg), Ptr{SpglibDB},
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    db = unsafe_load(ptr_cdb)
    return db
end

function get_spacegroup(lattice::Array{Float64, 2},
                        positions::Array{Float64, 2},
                        types::Array{Int64, 1},
                        symprec::Float64=1e-5)
    #
    num_atom = size(types)[1]

    db = spg_get_dataset(lattice, positions, types, num_atom, symprec)
    s = collect(Char.(db.international_symbol))
    return String(s), Base.convert(Int64, db.spacegroup_number)
end

function get_symmetry()
end

function get_smmetry_database()
end
