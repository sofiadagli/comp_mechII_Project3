"Main script for dynamic beam analysis"

import numpy as np
from assembly import assemble_global_matrices, apply_boundary_conditions
from static_solver import tip_load_vector, static_solve, analytical_solution

# Material & Geometry
E = 210e9
rho = 7850
L = 2.0
R = 0.02

A = np.pi * R**2
I = np.pi * R**4 /4.0

P = 1000.0
n_elem = 8

#--------------------------------------------------
# PART A: Static Analysis
#--------------------------------------------------
K, M = assemble_global_matrices(n_elem, E, I, rho, A, L)
f = tip_load_vector(n_elem, P)

fixed_dofs = [0, 1]

K_r, M_r, f_r, free = apply_boundary_conditions(K, M, f, fixed_dofs)

u = static_solve(K_r, f_r)

wL_num = u[-2]
wL_ana, _ = analytical_solution(P, L, E, I)

print(f"Static Tip Displacement FEM: {wL_num:.6e}")
print(f"Static Tip Displacement ANA: {wL_ana:.6e}")


