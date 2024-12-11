import os
import numpy as np
import sys
from mpi4py import MPI; rank = MPI.COMM_WORLD.Get_rank() # type: ignore
# from information import determine_next_filename
# from concurrent.futures import ProcessPoolExecutor

B_range = [.08,.08]


B_vec = np.linspace(*B_range,1)

def func(b):
    chain_length = (int(rank)+2)*2
    os.system(f'python julia_call.py {chain_length} {b}')

if __name__=="__main__":
    import time
    t0 = time.perf_counter()
    for i in B_vec:
        func(i)
    # with ProcessPoolExecutor(B_vec.size) as exe:
        # [_ for _ in exe.map(func, B_vec)]
    t1 = time.perf_counter()
    print(f'\n{np.round(t1-t0,3)} s')
