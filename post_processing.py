import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def animate_beam(
        u_hist,
        free_dofs,
        n_elem,
        L,
        time,
        scale=1.0,
        filename="beam.gif",
        fps0=30,
        max_frames=200,
        show=True
        ):
    
    n_steps, n_dof_free = u_hist.shape
    n_nodes = n_elem+1
    n_dof = 2*n_nodes
    
    step = max(1, n_steps // max_frames)
    u_hist = u_hist[::step]
    time = time[::step]
    
    #reconstract full dofs
    u_full = np.zeros((len(u_hist),n_dof))
    u_full[:, free_dofs] = u_hist
    
    w_hist = scale * u_full[:,0::2]
    x = np.linspace(0,L,n_nodes)
    
    fig, ax = plt.subplots()
    
    line, = ax.plot(x, w_hist[0], lw=2)
    
    #set axis limits
    w_max = np.max(np.abs(w_hist))
    ax.set_xlim(0,L)
    ax.set_ylim(-1.2*w_max, 1.2*w_max)
    
    time_text = ax.text(0.02,
                        0.90,
                        "",
                        transform=ax.transAxes,
                        fontsize=12,
                        bbox = dict(facecolor="green",alpha=0.8))
    
    def update(i):
        line.set_ydata(w_hist[i])
        time_text.set_text(f"t = {time[i]:.3f} s")
        
        return line, time_text
    
    ani = FuncAnimation(
        fig,
        update,
        frames=len(w_hist),
        interval = 1000/fps0,
        blit=True
        )
    
    ani.save(filename, writer=PillowWriter(fps=fps0))
    
    if show:
        plt.show()
    else:
        plt.close(fig)
        
    print(f"Animation saved as {filename}")
    
    
    
    


