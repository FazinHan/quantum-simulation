import sys
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

chain_length = int(sys.argv[1])

def layer_def(plan):
    if type(plan)==int:
        return '-'.join(['d']*plan)
    return '-'.join(['s' if j==1 else 'd' for j in plan])

for i in range(b_plot_count):
    
    b = B_vec[i+b_offset]
    
    name = determine_next_filename(f"julia_result_errors_{chain_length}_{b}_" ,'txt','data',exists=True)
    with open(name,'r') as file:
        # chain_length = eval(file.readline())
        layer_plan = eval(file.readline())
        Jii = eval(file.readline())
        fields = ''.join(file.readlines()).replace('  ',',').replace('\n',',')
        errors = 1 + np.array(eval(fields))
        # print(errors.shape);exit()
    x_axis = range(len(Jii))

    layer_list = np.array(range(len(layer_plan)))

    indices = np.argsort(errors[-1,:])[::-1]

    for idx, jii in enumerate(Jii):
        y_axis = errors[idx,:][indices]
        x_axis = layer_list[indices]
        axs.semilogy(layer_list, y_axis, label=f'$J_{{||}}={jii}$')

axs.set_xticks(layer_list, [layer_def(i) for i in layer_plan])
    
axs.grid()
axs.legend()
axs.set_xlabel('layers')
axs.set_ylabel('Abs. Error')
fig.suptitle(f'J = 1, $\\Omega={5*chain_length}$, B={B_vec[0]}, qubits={chain_length}')
name = determine_next_filename(f'jii_vs_err_{chain_length}_','png','outputs')
plt.savefig(name)
plt.clear()

### for the gradient calculation

fig, axs = plt.subplots(layout='constrained')

B_vec = np.linspace(*B_range, 1)

chain_length = int(sys.argv[1])

def layer_def(plan):
    if type(plan)==int:
        return '-'.join(['d']*plan)
    return '-'.join(['s' if j==1 else 'd' for j in plan])

for i in range(b_plot_count):
    
    b = B_vec[i+b_offset]
    
    name = determine_next_filename(f"julia_result_grads_{chain_length}_{b}_" ,'txt','data',exists=True)
    with open(name,'r') as file:
        # chain_length = eval(file.readline())
        layer_plan = eval(file.readline())
        Jii = eval(file.readline())
        grads = 1 + np.array(eval(file.readline()))
    x_axis = range(len(Jii))

    layer_list = np.array(range(len(layer_plan)))

    indices = np.argsort(grads[-1,:])[::-1]

    for idx, jii in enumerate(Jii):
        y_axis = grads[idx,:][indices]
        x_axis = layer_list[indices]
        axs.semilogy(layer_list, y_axis, label=f'$J_{{||}}={jii}$')

axs.set_xticks(layer_list, [layer_def(i) for i in layer_plan])
    
axs.grid()
axs.legend()
axs.set_xlabel('layers')
axs.set_ylabel('Gradients')
fig.suptitle(f'J = 1, $\\Omega={5*chain_length}$, B={B_vec[0]}, qubits={chain_length}')
name = determine_next_filename(f'jii_vs_grads_{chain_length}_','png','outputs')
plt.savefig(name)