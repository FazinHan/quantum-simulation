import os
import numpy as np
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor
from physics import B_range

B_vec = np.linspace(*B_range,2)

def func(b):
    os.system(f'python julia_call.py {b}')

if __name__=="__main__":
    # for i in B_vec:
        # func(i)
    with ProcessPoolExecutor(B_vec.size) as exe:
        [_ for _ in exe.map(func, B_vec)]