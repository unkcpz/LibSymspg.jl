export get_spacegroup, get_symmetry

mutable struct SpglibDB
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11, UInt8}
    hall_symbol::NTuple{17, UInt8}
    choice::NTuple{6, UInt8}
    transformation_matrix::NTuple{9, Float64}
    origin_shift::NTuple{3, Float64}
    n_operations::Cint
    rotations::Ptr{NTuple{9, Cint}}
    translations::Ptr{NTuple{3, Float64}}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{Tuple{7, UInt8}}
    equivalent_atoms::Ptr{Cint}
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{9, Float64}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3, Float64}}
    std_rotation_matrix::Ptr{NTuple{9, Float64}}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6, UInt8}
end

function spg_get_dataset(lattice::Array{Float64, 2},
                         positions::Array{Float64, 2},
                         types::Array{Int64, 1},
                         num_atom::Int64,
                         symprec::Float64=1e-5)

    @assert size(positions)[2] == size(types)[1] == num_atom

    types = Base.cconvert(Array{Int32, 1}, types)
    num_atom = Base.cconvert(Int32, num_atom)

    db =
    ccall( (:spg_get_dataset, libsymspg), Ref{SpglibDB},
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

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
    # BUG not output
    # println(collect(Char.(db.pointgroup_symbol)))
    return String(s), Base.convert(Int64, db.spacegroup_number)
end

"""
rotations output as 3x3xNop array
translations output as 3xNop array
"""
function get_symmetry(lattice::Array{Float64, 2},
                      positions::Array{Float64, 2},
                      types::Array{Int64, 1},
                      symprec::Float64=1e-5)
    #
    num_atom = size(types)[1]

    db = spg_get_dataset(lattice, positions, types, num_atom, symprec)
    nop = db.n_operations
    r = unsafe_wrap(Array{NTuple{9, Cint}}, db.rotations, nop)
    t = unsafe_wrap(Array{NTuple{3, Float64}}, db.translations, nop)
    r_ = Array{Int64, 3}(undef, 3, 3, nop)
    t_ = Array{Float64, 2}(undef, 3, nop)
    for i in 1:nop
        r_[:, :, i] = reshape([Base.convert(Int64, e) for e in r[i]], 3, 3)
        t_[:, i] = [e for e in t[i]]
    end

    n_atoms = db.n_atoms
    eq_atoms = unsafe_wrap(Array{Cint}, db.equivalent_atoms, n_atoms)
    return r_, t_, Base.convert(Array{Int64, 1}, eq_atoms)
end
