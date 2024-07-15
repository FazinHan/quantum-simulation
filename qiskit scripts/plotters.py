import matplotlib.pyplot as plt
import numpy as np

def classical_plotter(ax, B_arr, omega, energies, J, JII, plot_singlets=True):
    handles = []
    num_states = energies.shape[1]
    # print(num_states)
    for i in range(num_states):
        if i==0 and not plot_singlets:
            continue
        handles.append(ax.plot(B_arr/omega, energies[:,i],color='grey',ms=1.2,label=f'exact $|{i}\\rangle$')[0])
    return handles

def qiskit_plotter(ax, B_arr, energies, omega, J, JII, plot_singlets=True):
    # qiskit result plotter
    handles = []
    num_states = energies.shape[1]
    for i in range(num_states):
        if i==0 and not plot_singlets:
            continue
        handles.append(ax.plot(B_arr/omega, energies[:,i],'g.',label=f'VQE $|{i}\\rangle$')[0])
    return handles
    
def qiskit_cost_plotter(B_arr, costs, omega, J, JII):
    # costs plotter
    fig2, axs = plt.subplots(costs.shape[1],1) # number of layers is determined
    for idx, ax in enumerate(axs):        
        ax.plot(B_arr/omega, -costs[:,idx],'.')
        ax.set_xlabel('$B/\\Omega$')
        ax.set_ylabel('costs')# of $|{idx}\\rangle$')
    fig2.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    fig2.set_figheight(3.5*costs.shape[1])
    fig2.tight_layout()