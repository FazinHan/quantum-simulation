import os
import numpy as np
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor

B_range = [.08,.08]


B_vec = np.linspace(*B_range,1)

def func(b):
    chain_length = int(sys.argv[1])
    os.system(f'python julia_call.py {chain_length} {b}')

if __name__=="__main__":
    # for i in B_vec:
        # func(i)
    with ProcessPoolExecutor(B_vec.size) as exe:
        [_ for _ in exe.map(func, B_vec)]

