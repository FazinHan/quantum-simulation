import os
import matplotlib.pyplot as plt
import numpy as np

'''$(i) A$(i_A)_length$(chain_length)_depth$(depth) $(timeevo) $(en)'''

with open('julia_result.out','r') as file:
    data = file.readlines()
    data = [i.split() for i in data]

