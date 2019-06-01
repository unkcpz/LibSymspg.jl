export find_primitive, refine_cell, standardize_cell

"""
spg_standardize_cell standardize cell

lattice are represented as row vector here
means each colume is a basis of lattice
positions are represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function spg_standardize_cell(lattice::Array{Float64, 2},
                              positions::Array{Float64, 2},
                              types::Array{Int64, 1},
                              num_atom::Int64,
                              to_primitive::Int64,
                              no_idealize::Int64,
                              symprec::Float64=1e-5)

    @assert size(positions)[1] == size(types)[1] == num_atom

    # transpose to colomn vector used by spglib
    lattice_ = copy(lattice)
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
           lattice_, positions_, types_, num_atom, to_primitive, no_idealize, symprec )

    positions = positions_[:, 1:num_primitive_atom]
    types = types_[1:num_primitive_atom]

    # transpose back to row vectors
    positions = Array{Float64, 2}(positions')

    return lattice_, positions, Base.cconvert(Array{Int64, 1}, types), Base.cconvert(Int64,num_primitive_atom)
end

function standardize_cell(lattice::Array{Float64, 2},
                              positions::Array{Float64, 2},
                              types::Array{Int64, 1},
                              to_primitive::Bool,
                              no_idealize::Bool,
                              symprec::Float64=1e-5)

    to_primitive ? to_primitive_ = 1 : to_primitive_ = 0
    no_idealize ? no_idealize_ = 1 : no_idealize_ = 0

    num_atom = size(types)[1]

    return spg_standardize_cell(lattice, positions, types, num_atom, to_primitive_, no_idealize_, symprec)
end

"""
find_primitive reduce crystal to its primitive

lattice are represented as row vector here
means each colume is a basis of lattice
while positions are represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function find_primitive(lattice::Array{Float64, 2},
                            positions::Array{Float64, 2},
                            types::Array{Int64, 1},
                            symprec::Float64=1e-5)

    return standardize_cell(lattice, positions, types, true, false, symprec)
end

"""
refine_cell refine crystal to its convetional

lattice are represented as row vector here
means each colume is a basis of lattice
while positions are represented as row vector
means each row is a coordinate (x, y, z) of an atom
"""
function refine_cell(lattice::Array{Float64, 2},
                         positions::Array{Float64, 2},
                         types::Array{Int64, 1},
                         symprec::Float64=1e-5)

    return standardize_cell(lattice, positions, types, false, false, symprec)
end
