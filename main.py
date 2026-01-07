"Main script for dynamic beam analysis"

import numpy as np

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

