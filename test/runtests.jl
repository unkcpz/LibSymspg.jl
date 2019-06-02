using Test
using LibSymspg

latt = [4.0 0.0 0.0;
        0.0 4.0 0.0;
        0.0 0.0 4.0]
positions = [0.0 0.5; 0.0 0.5; 0.0 0.5]
types = [1, 1]
new_latt, new_positions, new_types = find_primitive(latt, positions, types, 1e-5)
# test arguments not modified
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.5; 0.0 0.5; 0.0 0.5]
@test types == [1, 1]

@test new_latt ≈ [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
@test new_positions ≈ [0.0, 0.0, 0.0]
@test new_types == [1]

latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
positions = Array{Float64, 2}([0.0 0.0 0.0]')
types = [1]
latt, positions, types = refine_cell(latt, positions, types, 1e-5)
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.5; 0.0 0.5; 0.0 0.5]
@test types == [1, 1]

# test determine the row and column type of latt and pos
latt = [4.0 0.0 0.0; 2.0 3.4641 0.0; 0.0 0.0 12.0]
positions = [0.0 1/3; 0.0 1/3; 0.0 1/3]
types = [1, 1]
num_atom = 2
s, n = get_spacegroup(latt, positions, types, 1e-3)
@test startswith(s, "P-3m1")
@test n == 164

r, t, eq_atoms = get_symmetry(latt, positions, types, 1e-3)
@test size(r)[3] == size(t)[2] == 12
@test eq_atoms == [0, 0]

latt = [4.0 20.0 0.0; 0.0 2.0 0.0; 0.0 0.0 12.0]
@test niggli_reduce(latt, 1e-3) ≈ [0.0 -2.0 0.0; 4.0 0.0 0.0; 0.0 0.0 12.0]
@test delaunay_reduce(latt, 1e-3) ≈ [0.0 2.0 0.0; -4.0 -0.0 0.0; -0.0 -0.0 12.0]

# ir_mesh test
latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
positions = Array{Float64, 2}([0.0 0.0 0.0]')
types = [1]
na = 1
meshk = [4, 4, 4]
is_shift = [0, 0, 0]
nir, kgrid, mapping = ir_reciprocal_mesh(latt, positions, types, na, meshk, is_shift, 1, 1e-5)
@test nir == 8
