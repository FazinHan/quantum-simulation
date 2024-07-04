import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
from optimiser import optimiser_main
from physics import Ω as omega
from physics import J, JII, B_range, num_rungs, num_layers, num_qubits
from information import determine_next_filename
import matplotlib.pyplot as plt
from plotters import qiskit_cost_plotter, qiskit_plotter, classical_plotter

if __name__=="__main__":

    B_size = 25

    B_arr = np.linspace(*B_range,B_size)*omega
    f_energies = []
    costs = []
    penalties = []
    ls = []

    with ProcessPoolExecutor(14) as exe:
        mapper = exe.map(optimiser_main, B_arr)
    for B, f_energy, cost, layer_step, penalty in mapper:
        f_energies.append(f_energy)
        costs.append(cost)
        ls.append(layer_step)
        penalties.append(penalty)
    # B_arr, singlets = zip(*sorted(singlets.items()))
    # B_arr, triplets = zip(*sorted(triplets.items()))
    f_energies = np.array(f_energies)
    penalties = np.array(penalties)
    ls = np.array(ls)
    costs = np.array(costs)

    
    filename = determine_next_filename('dimer','npz','data')
    with open(filename, 'wb') as file:
        np.savez(file, singlets=singlets, triplets=triplets, B_arr=B_arr, costs=costs, layer_step=layer_step, penalties=penalties)
        print('data saved in',filename)
    
    print(' ________ \n\n COMPLETE \n ________\n')
