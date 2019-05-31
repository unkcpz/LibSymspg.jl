using Test
using LibSymspg

latt = [4.0 0.0 0.0;
        0.0 4.0 0.0;
        0.0 0.0 4.0]
positions = [0.0 0.0 0.0; 0.5 0.5 0.5]
types = [1, 1]
num_atom = 2
new_latt, new_positions, new_types, new_num_atom = spg_find_primitive(latt, positions, types, num_atom, 1e-5)
# test arguments not modified
@test num_atom == 2
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.0 0.0; 0.5 0.5 0.5]
@test types == [1, 1]

@test new_num_atom == 1
@test new_latt ≈ [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
@test new_positions ≈ [0.0 0.0 0.0]
@test new_types == [1]

latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
positions = [0.0 0.0 0.0]
types = [1]
num_atom = 1
latt, positions, types, num_atom = spg_refine_cell(latt, positions, types, num_atom, 1e-5)
@test num_atom == 2
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.0 0.0; 0.5 0.5 0.5]
@test types == [1, 1]

latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
positions = [0.0 0.0 0.0]
types = [1]
num_atom = 1
