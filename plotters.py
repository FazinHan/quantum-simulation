import matplotlib.pyplot as plt
import numpy as np
from physics import B_range

def classical_plotter(ax, num_qubits, B_arr, omega, energies, J, JII, plot_singlets=True):
    for i in range(2**num_qubits):
        if i==0 and not plot_singlets:
            continue
        ax.plot(B_arr/omega, energies[:,i],'b--',ms=1.2,label='exact')
    # ax.set_xlim(*B_range)

def qiskit_plotter(ax, B_arr, singlets, triplets, omega, J, JII, plot_singlets=True):
    # qiskit result plotter
    if plot_singlets:
        ax.plot(B_arr/omega, singlets,'.',label='VQE')
    for i in range(3):
        ax.plot(B_arr/omega, triplets[:,i],'.',label='VQE')
    # ax.set_xlim(*B_range)
    
def qiskit_cost_plotter(ls, costs, omega, J, JII):
    # costs plotter
    fig2, axs = plt.subplots() # number of layers is determined
    axs.plot(ls, costs,'.')
    axs.set_xlabel('state label')
    axs.set_ylabel('costs')
    fig2.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    fig2.tight_layout()