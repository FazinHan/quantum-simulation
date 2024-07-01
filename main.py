import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename

B_arr = np.linspace(0,10,8)

def function(B):
    os.system(f'python optimiser.py {B}')

if __name__=="__main__":
    with ProcessPoolExecutor(8) as exe:
        mapper = exe.map(function , B_arr)
    
    for i in mapper:
        print(i)