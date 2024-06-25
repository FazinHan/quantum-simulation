from qutip import *
from physics import qutip_ladder as hamiltonian
import numpy as np

chain_length = 1
num_layers = 1
num_qubits = 2*chain_length

time = np.linspace(0,2*np.pi)

state = basis(num_qubits, 0)

sz = tensor(sigmaz(), sigmaz())

result = sesolve([[hamiltonian]], state, time, e_ops=[sz])