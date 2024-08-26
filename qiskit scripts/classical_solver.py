from qutip import *
from physics import qutip_ladder_hamiltonian as hamiltonian
from physics import floquetor
import numpy as np
# from physics import Ω as omega
from physics import J, JII, B_range, num_rungs, num_layers, num_qubits, T
from plotters import classical_plotter 
import matplotlib.pyplot as plt
from information import determine_next_filename

omega = 2.5

num_rungs = 2
num_layers = 2
num_qubits = 2*num_rungs

num_periods = 1
num_time_steps = 100
T = num_periods * 2*np.pi/omega
dt = T / num_time_steps

state = [basis(2,0)]*num_qubits
state[0] = basis(2,1)
state = tensor(state)

B_arr = np.linspace(0,.4,500)*omega

energies = []

for B in B_arr:

    ham_list = hamiltonian(num_rungs, B, omega, J, JII)
    f_basis = FloquetBasis(ham_list, T)
    
    # f_basis = floquetor(hamiltonian, T, B=B, omega=omega, rungs=num_rungs)
    f_energies = f_basis.e_quasi

    energies.append([*f_energies])

energies = np.array(energies)

# fname = determine_next_filename('qutip_data','npz','data')
# with open(fname,'wb') as file:
#     np.savez(file, B_arr=B_arr, energies=energies)
#     print(f'data saved to {fname}')

fig, ax = plt.subplots()
classical_plotter(ax, B_arr, omega, energies, J, JII)
fig.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
# plt.savefig(determine_next_filename('qutip',folder='outputs'))
plt.show()