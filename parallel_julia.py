import os
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor

Jii = np.arange(0,1,.1)

def func(jii):
    os.system(f'python julia_call.py {jii}')

if __name__=="__main__":
    with ProcessPoolExecutor(len(Jii)) as exe:
        [0 for _ in exe.map(func, Jii)]