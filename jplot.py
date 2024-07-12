import os
import matplotlib.pyplot as plt
import numpy as np
from physics import B_range
from julia_call import Jii
from information import determine_next_filename

b_offset = 0
b_plot_count = 1

layer_slice = '[:]'

fig, axs = plt.subplots(layout='constrained')


for i in range(b_plot_count):
    
    b = np.round(B_range[i+b_offset],1)

    name = determine_next_filename(f'julia_result_errors_0.0_','txt','data',exists=True)
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
