from qutip import *
from physics import qutip_ladder_hamiltonian as hamiltonian
from physics import floquetor
from information import determine_next_filename as filename
import numpy as np
import matplotlib.pyplot as plt
from physics import Ω as omega
from physics import J, JII

num_rungs = 1
num_layers = 1
num_qubits = 2*num_rungs

T = 2*np.pi/omega

state = [basis(2,0)]*num_qubits
state[0] = basis(2,1)
state = tensor(state)

B_arr = np.linspace(0,.1,1000)*omega

energies = []

for B in B_arr:

    ham_list = hamiltonian(num_rungs, B, omega, J, JII)
    f_basis = FloquetBasis(ham_list, T)
    
    # f_basis = floquetor(hamiltonian, T, B=B, omega=omega, rungs=num_rungs)
    f_energies = f_basis.e_quasi

    energies.append([*f_energies])

energies = np.array(energies)

for i in range(2**num_qubits):
    plt.plot(B_arr/omega, energies[:,i],'b.',ms=1.2)
plt.xlabel('$B/\\Omega$')
plt.ylabel('$\\varepsilon$')
plt.title(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
plt.savefig(filename())
# plt.show()