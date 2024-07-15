import os
import matplotlib.pyplot as plt
import numpy as np
from parallel_julia import B_range
from julia_call import Jii
from information import determine_next_filename

b_offset = 0
b_plot_count = 1

layer_slice = '[:]'

fig, axs = plt.subplots(layout='constrained')

B_vec = np.linspace(*B_range, 1)

def layer_def(plan):
    if type(plan)==int:
        return '-'.join(['d']*plan)
    return '-'.join(['s' if j==1 else 'd' for j in plan])

for i in range(b_plot_count):
    
    b = B_vec[i+b_offset]

    name = determine_next_filename(f'julia_result_errors_{b}_','txt','data',exists=True)
    with open(name,'r') as file:
        chain_length = eval(file.readline())
        layer_plan = eval(file.readline())
        Jii = eval(file.readline())
        errors = 1 + np.array(eval(file.readline()))
        print(errors)
    x_axis = range(len(Jii))

    for idx, layers in enumerate(layer_plan):
        axs.semilogy(eval('Jii'+layer_slice), eval('errors[:,idx]'+layer_slice), label=f'{layer_def(layers)}')
    
axs.set_xticks(x_axis, np.array(Jii))
axs.grid()
axs.legend()
axs.set_xlabel('$J||$')
axs.set_ylabel('Abs. Error')
fig.suptitle(f'J = 1, $\\Omega={5*chain_length}$, B={B_vec[0]}, qubits={chain_length}')
name = determine_next_filename(f'jii_vs_err_{chain_length}_','png','outputs')
plt.savefig(name)
# plt.show()