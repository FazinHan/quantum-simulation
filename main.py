import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from optimiser import optimiser_main
from information import determine_next_filename
from physics import Ω as omega

B_size = 50

B_arr = np.linspace(0,10,B_size)*omega
singlets = []
triplets = []

if __name__=="__main__":
    with ProcessPoolExecutor(B_size) as exe:
        mapper = exe.map(optimiser_main, B_arr)
    for singlet, triplet in mapper:
        singlets.append(singlet)
        triplets.append(triplet)
    plt.plot(B_arr, singlets,'.')
    triplets = np.array(triplets)
    for i in range(3):
        plt.plot(B_arr, triplets[:,i],'.')
    plt.xlabel('$B/\\Omega$')
    plt.ylabel('$\\epsilon$')
    plt.savefig(determine_next_filename())

    with open(determine_next_filename('data','npz'), 'wb') as file:
        np.savez(file, singlets=singlets, triplets=triplets, B_arr=B_arr, B_size=B_size)

    print(' ________ \n\n COMPLETE \n ________\n')