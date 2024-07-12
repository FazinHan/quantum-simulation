import os
import matplotlib.pyplot as plt
import numpy as np
from parallel_julia import B
from julia_call import Jii
from information import determine_next_filename

b_offset = 1
b_plot_count = 4

layer_slice = '[:]'

fig, axs = plt.subplots(layout='constrained')


for i in range(b_plot_count):
    
    b = np.round(B[i+b_offset],1)
    
    name = determine_next_filename(f'julia_result_errors0.4','txt','data//round1',exists=True)
    with open(name,'r') as file:
        layer_plan = eval(file.readline())
        B = eval(file.readline())
        errors = 1 + np.array(eval(file.readline()))
    x_axis = range(len(layer_plan))
    
    axs.semilogy(eval('x_axis'+layer_slice), eval('errors[i]'+layer_slice), label=f'B={b}')
    
axs.set_xticks(x_axis, layer_plan)
axs.grid()
axs.legend()
fig.suptitle('J = 1')
plt.show()
