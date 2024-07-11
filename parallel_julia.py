import os
import numpy as np
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor

Jii = np.arange(0,1,.1)

def func(jii):
    os.system(f'python julia_call.py {jii}')

if __name__=="__main__":
    for i in Jii:
        func(i)