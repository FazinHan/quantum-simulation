import matplotlib.pyplot as plt
import numpy as np
from information import determine_next_filename
from physics import Ω as omega
from physics import J, JII, B_range, num_rungs, num_layers, num_qubits
from plotters import *

if __name__=="__main__":
    fig, ax = plt.subplots()
    
    with open(determine_next_filename('qutip_data','npz','data',exists=True),'rb') as file:
        data = np.load(file)
        B_arr = data['B_arr']
        energies = data['energies']
        # print(B_arr[-1])
    classical_plotter(ax, num_qubits, B_arr, omega, energies, J, JII, plot_singlets=False)

    with open(determine_next_filename('dimer','npz','data',exists=True),'rb') as file:
        data = np.load(file)
        B_arr = data['B_arr']
        singlets = data['singlets']
        triplets = data['triplets']
        ls = data['layer_step']
        # print(B_arr[-1])
        costs = data['costs']
        penalties = data['penalties']
    qiskit_plotter(ax, B_arr, singlets, triplets, omega, J, JII)
    
    ax.set_xlabel('$B/\\Omega$')
    ax.set_ylabel('$\\epsilon$')
    # ax.set_xlim(0,.1)
    fig.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    plt.savefig(determine_next_filename(folder='outputs'))

    qiskit_cost_plotter(ls, costs, omega, J, JII)
    # qiskit_penalty_plotter(ls, penalties, omega, J, JII)
    plt.savefig(determine_next_filename(folder='outputs'))