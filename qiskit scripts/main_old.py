import matplotlib.pyplot as plt
import numpy as np
from information import determine_next_filename
from physics import Ω as omega
from physics import J, JII, B_range, num_rungs, num_layers, num_qubits
from plotters import *
import sys



if __name__=="__main__1":
    with open(determine_next_filename('dimer','npz','data',exists=True),'rb') as file:
        data = np.load(file)
        B_arr = data['B_arr']
        # singlets = data['singlets']
        # triplets = data['triplets']
        ls = data['layer_step']
        qenergies = data['f_energies']
        costs = data['costs']
        penalties = data['penalties']
    handles = qiskit_plotter(ax, B_arr, qenergies, omega, J, JII)

if __name__=="__main__":
    fig, ax = plt.subplots()

    with open(determine_next_filename('qutip_data','npz','data',exists=True),'rb') as file:
        data = np.load(file)
        B_arr = data['B_arr']
        qenergies = data['energies']
    classical_plotter(ax, B_arr, omega, qenergies, J, JII, plot_singlets=True)
    ax.set_xlabel('$B/\\Omega$')
    ax.set_ylabel('$\\epsilon$')
    # ax.set_xlim(0,.1)
    # ax.legend(handles=handles)
    fig.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    fig.tight_layout()
    
    name = determine_next_filename(fname='q_energies',folder='outputs')
    plt.savefig(name)
    print(name, 'created')


if __name__=="__main__1":

    qiskit_cost_plotter(B_arr, costs, omega, J, JII)
    # qiskit_cost_plotter(ls, penalties, omega, J, JII)
    name = determine_next_filename(fname='costs',folder='outputs')
    plt.savefig(name)
    print(name, 'created')