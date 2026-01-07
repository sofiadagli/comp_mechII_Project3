"Newmark time integration method"

import numpy as np

def newmark(M, K, f_func, u0, v0, dt, t_end, alpha=0.25, beta=0.5):
    n = len(u0)
    n_steps = int(t_end/dt)
    t = np.linspace(0, t_end, n_steps)
    
    #init
    u = np.zeros((n_steps,n))
    v = np.zeros((n_steps,n))
    a = np.zeros((n_steps,n))
    
    u[0] = u0
    v[0] = v0
    a[0] = np.linalg.solve(M, f_func(0)-K@u0)
    
    K_eff = K + M/(alpha*dt**2)
    K_eff_inv = np.linalg.inv(K_eff)
    
    for i in range(1,n_steps):
        f_eff = (
            f_func(t[i]) + 
            M @ (
                u[i-1]/ (alpha*dt**2) +
                v[i-1]/ (alpha*dt) +
                (1/(2*alpha)-1)*a[i-1]
                ))
        
        u[i] = K_eff_inv @ f_eff
        a[i] = (
            (u[i] - u[i-1]) / (alpha*dt**2)
            - v[i-1] / (alpha*dt)
            - (1/(2*alpha)-1) * a[i-1]
            )
        v[i] = v[i-1] + dt*((1-beta)*a[i-1]+beta*a[i])
        
    return t, u, v, a

        
    
    
    

