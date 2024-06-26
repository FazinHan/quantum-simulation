from qutip import *
from physics import qutip_ladder_hamiltonian as hamiltonian
from information import determine_next_filename as filename
import numpy as np
import matplotlib.pyplot as plt

num_rungs = 1
num_layers = 1
num_qubits = 2*num_rungs

omega = .5

T = 2*np.pi/omega

state = [basis(2,0)]*num_qubits
state[0] = basis(2,1)
state = tensor(state)

B_arr = np.linspace(0,5,500)

energies = []

for B in B_arr:
    
    ham_list = hamiltonian(num_rungs, J=2)

    # print(ham_list[0])
    
    f_basis = FloquetBasis(ham_list, T, args={'B':B*omega, 'omega':omega})
    f_energies = f_basis.e_quasi

    energies.append([*f_energies])

energies = np.array(energies)

for i in range(2**num_qubits):
    plt.plot(B_arr, energies[:,i])
plt.xlabel('$B/\\Omega$')
plt.ylabel('$\\varepsilon$')
plt.savefig(filename())