import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from optimiser import optimiser_main
from information import determine_next_filename

B_arr = np.linspace(0,10,8)
singlets = []
triplets = []

if __name__=="__main__":
    with ProcessPoolExecutor(8) as exe:
        mapper = exe.map(optimiser_main, B_arr)
    for singlet, triplet in mapper:
        singlets.append(singlet)
        triplets.append(triplet)
    plt.plot(B_arr, singlets,'.')
    triplets = np.array(triplets)
    for i in range(3):
        plt.plot(B_arr, triplets[:,i],'.')
    plt.savefig(determine_next_filename())