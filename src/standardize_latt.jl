export niggli_reduce, delaunay_reduce

function niggli_reduce(lattice::Array{Float64, 2},
                       symprec::Float64=1e-5)
    #
    lattice_ = copy(lattice)
    r = ccall((:spg_niggli_reduce, libsymspg), Int32,
              (Ptr{Float64}, Float64),
              lattice_, symprec)
    @assert r != 0
    return lattice_
end

function delaunay_reduce(lattice::Array{Float64, 2},
                         symprec::Float64=1e-5)
    #
    lattice_ = copy(lattice)
    r = ccall((:spg_delaunay_reduce, libsymspg), Int32,
              (Ptr{Float64}, Float64),
              lattice_, symprec)
    @assert r != 0
    return lattice_
end
