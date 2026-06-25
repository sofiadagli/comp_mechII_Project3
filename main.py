"""
Main script for dynamic beam analysis
"""

import numpy as np
import matplotlib.pyplot as plt

from assembly import assemble_global_matrices, apply_boundary_conditions
from static_solver import tip_load_vector, static_solve, analytical_solution
from newmark import newmark, freq_estimation
from post_processing import animate_beam

# ============================================================
# Material & Geometry
# ============================================================

E = 210e9
rho = 7850
L = 2.0
R = 0.02

A = np.pi * R**2
I = np.pi * R**4 / 4.0

P = 1000.0
n_elem = 26

fixed_dofs = [0, 1]

# ============================================================
# MESH CONVERGENCE STUDY
# ============================================================

mesh_list = [1, 2, 4, 8, 16, 32]
errors = []

wL_exact, theta_exact = analytical_solution(P, L, E, I)

for ne in mesh_list:

    K_tmp, M_tmp = assemble_global_matrices(
        ne, E, I, rho, A, L
    )

    f_tmp = tip_load_vector(ne, P)

    K_r_tmp, M_r_tmp, f_r_tmp, free_tmp = \
        apply_boundary_conditions(
            K_tmp,
            M_tmp,
            f_tmp,
            fixed_dofs
        )

    u_tmp = static_solve(K_r_tmp, f_r_tmp)

    wL_num_tmp = u_tmp[-2]

    err = abs(wL_num_tmp - wL_exact) / abs(wL_exact)

    errors.append(err)

plt.figure()

plt.loglog(mesh_list, errors, "o-")

plt.xlabel("Number of elements")
plt.ylabel("Relative error")
plt.title("Mesh convergence")

plt.grid(True)

plt.show()

# ============================================================
# PART A: STATIC ANALYSIS
# ============================================================

K, M = assemble_global_matrices(
    n_elem,
    E,
    I,
    rho,
    A,
    L
)

f = tip_load_vector(n_elem, P)

K_r, M_r, f_r, free = apply_boundary_conditions(
    K,
    M,
    f,
    fixed_dofs
)

u = static_solve(K_r, f_r)

wL_num = u[-2]
theta_num = u[-1]

wL_ana, theta_ana = analytical_solution(
    P,
    L,
    E,
    I
)

print("\n===== STATIC ANALYSIS =====")

print(f"Tip displacement FEM : {wL_num:.6e}")
print(f"Tip displacement ANA : {wL_ana:.6e}")

print(f"Tip rotation FEM     : {theta_num:.6e}")
print(f"Tip rotation ANA     : {theta_ana:.6e}")

# ============================================================
# PART B: FREE VIBRATION
# ============================================================

u0 = u.copy()
v0 = np.zeros_like(u0)

zero_force = lambda t: np.zeros_like(u0)

dt = 1e-4
t_end = 2.8

t, uh, vh, ah = newmark(
    M_r,
    K_r,
    zero_force,
    u0,
    v0,
    dt,
    t_end
)

omega_1 = freq_estimation(
    t,
    uh[:, -2],
    min_height=0.0,
    min_distance=0.08
)

omega_ana = (
    1.875**2
) * np.sqrt(
    E * I /
    (rho * A * L**4)
)

print("\n===== FREE VIBRATION =====")

print(f"ω1 FEM : {omega_1:.3f}")
print(f"ω1 ANA : {omega_ana:.3f}")

# ============================================================
# PART C: FORCED VIBRATION
# ============================================================

P0 = 1000.0

Omega = 0.95 * omega_1

def harmonic_force(t):

    f = np.zeros_like(u0)

    f[-2] = P0 * np.sin(Omega * t)

    return f

u0 = np.zeros_like(u0)
v0 = np.zeros_like(u0)

t, uh, vh, ah = newmark(
    M_r,
    K_r,
    harmonic_force,
    u0,
    v0,
    dt,
    t_end
)

# ============================================================
# TIP DISPLACEMENT
# ============================================================

plt.figure()

plt.plot(t, uh[:, -2])

plt.xlabel("Time (s)")
plt.ylabel("Tip displacement (m)")
plt.title("Free-end displacement")

plt.grid(True)

plt.show()

# ============================================================
# MIDSPAN DISPLACEMENT
# ============================================================

mid_node = n_elem // 2

mid_dof_global = 2 * mid_node

mid_index = np.where(
    free == mid_dof_global
)[0][0]

plt.figure()

plt.plot(
    t,
    uh[:, mid_index]
)

plt.xlabel("Time (s)")
plt.ylabel("Midspan displacement (m)")
plt.title("Midspan response")

plt.grid(True)

plt.show()

# ============================================================
# PHASE DIAGRAM
# ============================================================

plt.figure()

plt.plot(
    uh[:, -2],
    vh[:, -2]
)

plt.xlabel("Displacement (m)")
plt.ylabel("Velocity (m/s)")
plt.title("Phase diagram")

plt.grid(True)

plt.show()

# ============================================================
# ENERGIES
# ============================================================

kinetic = np.zeros(len(t))
strain = np.zeros(len(t))

for i in range(len(t)):

    kinetic[i] = (
        0.5 *
        vh[i] @ M_r @ vh[i]
    )

    strain[i] = (
        0.5 *
        uh[i] @ K_r @ uh[i]
    )

total = kinetic + strain

# ============================================================
# ENERGY PLOT
# ============================================================

plt.figure()

plt.plot(
    t,
    kinetic,
    label="Kinetic Energy"
)

plt.plot(
    t,
    strain,
    label="Strain Energy"
)

plt.plot(
    t,
    total,
    label="Total Energy"
)

plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.title("Energy evolution")

plt.legend()

plt.grid(True)

plt.show()

# ============================================================
# ANIMATION
# ============================================================

animate_beam(
    uh,
    free_dofs=free,
    n_elem=n_elem,
    L=L,
    time=t,
    scale=55,
    filename="semfebeam.gif",
    fps0=30,
    max_frames=150,
    show=True
)