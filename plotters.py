import matplotlib.pyplot as plt
import numpy as np
from physics import B_range

def classical_plotter(ax, B_arr, omega, energies, J, JII, plot_singlets=True):
    handles = []
    num_states = np.unique(energies[:,0]).size
    for i in range(num_states):
        if i==0 and not plot_singlets:
            continue
        handles.append(ax.plot(B_arr/omega, energies[:,np.where(energies[:,0]==i)],color='grey',ms=1.2,label=f'exact $|{i}\\rangle$')[0])
    return handles

def qiskit_plotter(ax, B_arr, energies, omega, J, JII, plot_singlets=True):
    # qiskit result plotter
    handles = []
    num_states = np.unique(energies[:,0]).size
    for i in range(num_states):
        if i==0 and not plot_singlets:
            continue
        handles.append(ax.plot(B_arr/omega, energies[:,np.where(energies[:,0]==i)],'.',label=f'VQE $|{i}\\rangle$')[0])
    return handles
    
def qiskit_cost_plotter(B_arr, costs, omega, J, JII):
    # costs plotter
    fig2, axs = plt.subplots(costs.shape[1],1) # number of layers is determined
    for idx, ax in enumerate(axs):        
        ax.plot(B_arr/omega, costs[:,idx],'.')
        ax.set_xlabel('$B/\\Omega$')
        ax.set_ylabel(f'costs of $|{idx}\\rangle$')
    fig2.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    fig2.set_figheight(7*2)
    fig2.tight_layout()