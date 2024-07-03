import matplotlib.pyplot as plt
import numpy as np

def classical_plotter(ax, num_qubits, B_arr, omega, energies, J, JII):
    for i in range(2**num_qubits):
        ax.plot(B_arr/omega, energies[:,i],'b.',ms=1.2)
    

def qiskit_plotter(ax, B_arr, singlets, triplets, omega, J, JII):
    # qiskit result plotter
    ax.plot(B_arr, singlets,'.')
    for i in range(3):
        ax.plot(B_arr, triplets[:,i],'.')
    
def qiskit_cost_plotter(B_arr, ls, costs, omega, J, JII):
    # costs plotter
    fig2, axs = plt.subplots(np.unique(ls[:,0]).size,1) # number of layers is determined
    for idx, n in enumerate(ls[:,0]):
        locs = np.where(ls[:,0]==n)
        try:
            axs[idx].plot(ls[locs][:,1], costs[locs],'.')
            axs[idx].set_xlabel(f'layer {n[0]} steps')
            axs[idx].set_ylabel('cost')
        except TypeError:
            axs.plot(ls[locs][:,1], costs[locs],'.')
            axs.set_xlabel('steps')
            axs.set_ylabel('costs')
    fig2.tight_layout()
    fig2.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    