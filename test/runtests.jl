using Test
using LibSymspg

latt = [4.0 0.0 0.0;
        0.0 4.0 0.0;
        0.0 0.0 4.0]
positions = [0.0 0.0 0.0; 0.5 0.5 0.5]
types = [1, 1]
new_latt, new_positions, new_types = find_primitive(latt, positions, types, 1e-5)
# test arguments not modified
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.0 0.0; 0.5 0.5 0.5]
@test types == [1, 1]

@test new_latt ≈ [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
@test new_positions ≈ [0.0 0.0 0.0]
@test new_types == [1]

latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
positions = [0.0 0.0 0.0]
types = [1]
latt, positions, types = refine_cell(latt, positions, types, 1e-5)
@test latt ≈ [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]
@test positions ≈ [0.0 0.0 0.0; 0.5 0.5 0.5]
@test types == [1, 1]

# latt = [-2.0 2.0 2.0; 2.0 -2.0 2.0; 2.0 2.0 -2.0]
# positions = [0.0 0.0 0.0]
# types = [1]
# num_atom = 1
# db = spg_get_dataset(latt, positions, types, num_atom, 1e-5)
