from qutip import *
from physics import qutip_ladder_hamiltonian as hamiltonian
from physics import floquetor
import numpy as np
from physics import Ω as omega
from physics import J, JII, B_range
from plotters import classical_plotter 
import matplotlib.pyplot as plt
from information import determine_next_filename

num_rungs = 1
num_layers = 1
num_qubits = 2*num_rungs

T = 2*np.pi/omega

state = [basis(2,0)]*num_qubits
state[0] = basis(2,1)
state = tensor(state)

B_arr = np.linspace(*B_range,1000)*omega

energies = []

for B in B_arr:

    ham_list = hamiltonian(num_rungs, B, omega, J, JII)
    f_basis = FloquetBasis(ham_list, T)
    
    # f_basis = floquetor(hamiltonian, T, B=B, omega=omega, rungs=num_rungs)
    f_energies = f_basis.e_quasi

    energies.append([*f_energies])

energies = np.array(energies)

with open() as file:
    np.savez(file, B_arr=B_arr, f_energies=f_energies)

classical_plotter(num_qubits, B_arr, omega, energies, J, JII)
plt.savefig(determine_next_filename())
# plt.show()