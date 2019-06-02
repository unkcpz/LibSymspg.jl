export ir_reciprocal_mesh

function ir_reciprocal_mesh(lattice::Array{Float64, 2},
                            positions::Array{Float64, 2},
                            types::Array{Int64, 1},
                            num_atom::Int64,
                            meshk::Array{Int64, 1},
                            is_shift::Array{Int64},
                            is_time_reversal=1,
                            symprec::Float64=1e-5)
    #
    @assert size(positions)[2] == size(types)[1] == num_atom

    cmeshk = Base.cconvert( Array{Int32,1}, meshk )
    cis_shift = Base.cconvert( Array{Int32,1}, is_shift )
    ctypes = Base.cconvert( Array{Int32,1}, types)
    num_atom = Base.cconvert( Int32, num_atom )
    is_t_rev = Base.cconvert( Int32, is_time_reversal )

    # Prepare for output
    Nkpts = prod(meshk)
    kgrid = zeros(Int32,3,Nkpts)
    mapping = zeros(Int32,Nkpts)

    num_ir =
    ccall((:spg_get_ir_reciprocal_mesh, libsymspg), Int32,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Int32, Ptr{Float64}, Ptr{Float64},
           Ptr{Int32}, Int32, Float64),
           kgrid, mapping, cmeshk, cis_shift, is_t_rev,
           lattice, positions, ctypes, num_atom, symprec)

    return Base.cconvert(Int64, num_ir),
           Base.cconvert(Array{Int64,2}, kgrid),
           Base.cconvert(Array{Int64,1}, mapping)
end
