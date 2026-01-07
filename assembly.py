"Assembly of the global matrices"

import numpy as np

def beam_element_matrices(E, I, rho, A, Le):
    ke = (E*I/Le**3)*np.array([
        [  12,    6*Le,   -12,    6*Le]
        [6*Le, 4*Le**2, -6*Le, 2*Le**2]
        [ -12,   -6*Le,    12,   -6*Le]
        [6*Le, 2*Le**2, -6*Le, 4*Le**2]])
    
    me = (rho*A*Le/420)*np.array([
        [   156,    22*Le,     54,   -13*Le]
        [ 22*Le,  4*Le**2,  13*Le, -3*Le**2]
        [    54,    13*Le,    156,   -22*Le]
        [-13*Le, -3*Le**2, -22*Le,  4*Le**2]])
    
    return ke, me


def assemble_global_matrices(n_elem, E, I, rho, A, L):
    n_nodes = n_elem + 1
    dof = 2*n_nodes
    Le = L/n_elem
    
    K = np.zeros((dof, dof))
    M = np.zeros((dof, dof))

    for e in range(n_elem):
        ke, me = beam_element_matrices(E, I, rho, A, Le)
        idx = np.array([2*e,2*e+1,2*e+2,2*e+3])
        
        for i in range(4):
            for j in range(4):
                K[idx[i], idx[j]] += ke[i,j]
                M[idx[i], idx[j]] += me[i,j]
                
    return K, M

                
        
        
    