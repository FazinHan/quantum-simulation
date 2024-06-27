import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename

B_arr = np.linspace(0,10,8)

with ProcessPoolExecutor(4) as exe:
    mapper = exe.map(lambda B: os.system(f'python optimiser.py {B}'), B_arr)

_ = [i for i in mapper]