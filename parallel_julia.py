import os
import sys
from information import determine_next_filename
from concurrent.futures import ProcessPoolExecutor

Jii = range(1,11)

def func(jii):
    os.system(f'python julia_call.py {int(jii)}')

if __name__=="__main__":
    with ProcessPoolExecutor(len(Jii)) as exe:
        [0 for _ in exe.map(func, Jii)]