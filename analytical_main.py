from qutip import *
from physics import qutip_ladder_hamiltonian as hamiltonian
import numpy as np

chain_length = 1
num_layers = 1
num_qubits = 2*chain_length

time = np.linspace(0,2*np.pi)

state = basis(num_qubits, 0)

sz = tensor(sigmaz(), sigmaz())

ham_list = hamiltonian(1)

print(ham_list)

result = sesolve(ham_list, state, time, e_ops=[sz])