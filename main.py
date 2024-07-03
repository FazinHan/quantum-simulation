import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
from optimiser import optimiser_main
from physics import Ω as omega
from physics import J, JII, B_range, num_rungs, num_layers, num_qubits
from information import determine_next_filename
import matplotlib.pyplot as plt
from plotters import qiskit_cost_plotter, qiskit_plotter, classical_plotter

B_size = 25

B_arr = np.linspace(*B_range,B_size)*omega
singlets = {}
triplets = {}
costs = []
ls = []

if __name__=="__main__":
    with ProcessPoolExecutor(10) as exe:
        mapper = exe.map(optimiser_main, B_arr)
    for B, singlet, triplet, cost, layer_step in mapper:
        singlets[B] = singlet
        triplets[B] = triplet
        costs.append(cost)
        ls.append(layer_step)
    B_arr, singlets = zip(*sorted(singlets.items()))
    B_arr, triplets = zip(*sorted(triplets.items()))
    singlets = np.array(singlets)
    triplets = np.array(triplets)
    ls = np.array(ls)
    costs = np.array(costs)

    
    filename = determine_next_filename('dimer','npz')
    with open(filename, 'wb') as file:
        np.savez(file, singlets=singlets, triplets=triplets, B_arr=B_arr, costs=costs, layer_step=layer_step)
        print('data saved in',filename)

    qiskit_plotter(B_arr, singlets, triplets, omega, J, JII)
    with open(determine_next_filename('qutip_data','npz',exists=True),'rb') as file:
        data = np.load(file)
        B_arr = data['B_arr']
        energies = data['energies']
    classical_plotter(num_qubits, B_arr, omega, energies, J, JII)
    plt.savefig(determine_next_filename())

    qiskit_cost_plotter(B_arr, ls, costs, omega, J, JII)
    plt.savefig(determine_next_filename())
    
    print(' ________ \n\n COMPLETE \n ________\n')