import os
import numpy as np
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor
from physics import B_range

def func(b):
    os.system(f'python julia_call.py {b}')

if __name__=="__main__":
    for i in B_range:
        func(i)