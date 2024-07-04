import matplotlib.pyplot as plt
import numpy as np
from physics import B_range

def classical_plotter(ax, num_qubits, B_arr, omega, energies, J, JII, plot_singlets=True):
    handles = []
    for i in range(2**num_qubits):
        if i==0 and not plot_singlets:
            continue
        handles.append(ax.plot(B_arr/omega, energies[:,i],color='grey',ms=1.2,label='exact'))
    return handles
    # ax.set_xlim(*B_range)

def qiskit_plotter(ax, B_arr, singlets, triplets, omega, J, JII, plot_singlets=True):
    # qiskit result plotter
    handles = []
    if plot_singlets:
        handles.append(ax.plot(B_arr/omega, singlets,'.',label='VQE'))
    for i in range(3):
        handles.append(ax.plot(B_arr/omega, triplets[:,i],'.',label='VQE'))
    # ax.set_xlim(*B_range)
    return handles
    
def qiskit_cost_plotter(B_arr, costs, omega, J, JII):
    # costs plotter
    fig2, axs = plt.subplots(costs.shape[1],1) # number of layers is determined
    for idx, ax in enumerate(axs):        
        ax.plot(B_arr/omega, costs[:,idx],'.')
        ax.set_xlabel('$B/\\Omega$')
        ax.set_ylabel('costs')
    fig2.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    fig2.set_figheight(7*2)
    fig2.tight_layout()