import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename

B_arr = np.linspace(0,10,8)

if __name__=="__main__":
    with ProcessPoolExecutor(4) as exe:
        mapper = exe.map(lambda B: os.system(f'python optimiser.py {B}'), B_arr)
    
    for i in mapper:
        print(i)