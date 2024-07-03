import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from optimiser import optimiser_main
from information import determine_next_filename
from physics import Ω as omega
from physics import J, JII

B_size = 100

B_arr = np.linspace(0,10,B_size)*omega
singlets = {}
triplets = {}
costs = []

if __name__=="__main__":
    with ProcessPoolExecutor(10) as exe:
        mapper = exe.map(optimiser_main, B_arr)
    for B, singlet, triplet, cost, layer_step in mapper:
        singlets[B] = singlet
        triplets[B] = triplet
        costs.append(cost)
    B_arr, singlets = zip(*sorted(singlets.items()))
    B_arr, triplets = zip(*sorted(triplets.items()))
    plt.plot(B_arr, singlets,'.')
    triplets = np.array(triplets)
    for i in range(3):
        plt.plot(B_arr, triplets[:,i],'.')
    plt.xlabel('$B/\\Omega$')
    plt.ylabel('$\\epsilon$')
    plt.savefig(determine_next_filename())

    with open(determine_next_filename('dimer','data','npz'), 'wb') as file:
        np.savez(file, singlets=singlets, triplets=triplets, B_arr=B_arr, costs=costs, layer_step=layer_step)

    ls, costs = zip(*sorted(costs.items()))
    fig, axs = plt.subplots(np.unique(ls[:,0]).size,1)
    fig.suptitle(f'$\\Omega={omega}$, $J={J}$, $J||={JII}$')
    for idx, n in enumerate(np.unique(ls[:,0])):
        locs = np.where(ls[:,0]==n)
        axs[idx].plot(ls[locs][:,1], costs[locs])
    fig.tight_layout()
    plt.savefig(determine_next_filename())
    
    
    

    print(' ________ \n\n COMPLETE \n ________\n')