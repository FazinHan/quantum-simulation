import os
import matplotlib.pyplot as plt
import numpy as np
from parallel_julia import JII
from information import determine_next_filename

jii_offset = 1
jii_plot_count = 4

layer_slice = '[:]'

fig, axs = plt.subplots(layout='constrained')


for i in range(jii_plot_count):
    
    jii = np.round(JII[i+jii_offset],1)
    
    name = determine_next_filename(f'julia_result_errors{jii}','txt','data',exists=True)
    with open(name,'r') as file:
        layer_plan = eval(file.readline())
        Jii = eval(file.readline())
        errors = 1 + np.array(eval(file.readline()))
    x_axis = range(len(layer_plan))
    
    axs.plot(eval('x_axis'+layer_slice), eval('errors[i]'+layer_slice), label=f'J||={jii}')
    
axs.set_xticks(x_axis, layer_plan)
axs.grid()
axs.legend()
fig.suptitle('J = 1')
plt.show()
