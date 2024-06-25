from qiskit.quantum_info import SparsePauliOp
import numpy as np
from qiskit.circuit.library import HamiltonianGate

Ω = 2.5
h_cut = 1

def hamiltonian_circular(t, A=2, J=1, omega=Ω):
    creator = ['I']*chain_length
    paulis = ['I','X','Y','Z']
    ham = [] # [('X',1.0)]
    for i in range(chain_length-1):
        for j in range(1,4):
            op = creator[:]
            op[i] = paulis[j]
            op[i+1] = paulis[j]
            ham.append([''.join(op), -J/4])
    for i in range(chain_length):
        op1, op2 = creator[:], creator[:]
        op1[i] = 'X'
        op2[i] = 'Y'
        ham.append([''.join(op1), A * np.cos(omega*t)])
        ham.append([''.join(op2), A * np.sin(omega*t)])
    ham = np.array(ham)
    # print(A * np.cos(Ω*t))
    return SparsePauliOp(ham[:,0], ham[:,1])

def hamiltonian_linear(t, A, Δ=1, omega=Ω):
    ham = SparsePauliOp(['Z','X'] , [-Δ/2, A/2*np.cos(omega*t)])
    # plt.plot(t, A*np.cos(Ω*t)/2,'.')
    return ham


def hamiltonian_ladder(t, num_states, entanglement_ratio, J=1, B=1, omega=2.5):
    J11 = entanglement_relator(J, entanglement_ratio)
    pauli = ['I','X','Y','Z']
    ham = []
    coeffs = []
    for i in range(int(num_states/2)):
        creator = ['I']*num_states
        creator[i] = pauli[1]
        coeff = B * np.cos(omega * t)
        ham.append(''.join(creator))
        coeffs.append(coeff)
        creator = ['I']*num_states
        creator[i+1] = pauli[1]
        coeff = B * np.cos(omega * t)
        ham.append(''.join(creator))
        coeffs.append(coeff)
        creator = ['I']*num_states
        creator[i] = pauli[2]
        coeff = B * np.sin(omega * t)
        ham.append(''.join(creator))
        coeffs.append(coeff)
        creator = ['I']*num_states
        creator[i+1] = pauli[2]
        coeff = B * np.sin(omega * t)
        ham.append(''.join(creator))
        coeffs.append(coeff)
        for j in range(1,4):
            creator = ['I']*num_states
            creator[2*i] = pauli[j]
            creator[2*i+1] = pauli[j]
            ham.append(''.join(creator))
            coeffs.append(J)
    for i in range(int(num_states/2)-1):
        for j in range(1,4):
            creator = ['I']*num_states
            creator[2*i] = pauli[j]
            creator[2*(i+1)] = pauli[j]
            ham.append(''.join(creator))
            coeffs.append(J11)
            creator = ['I']*num_states
            creator[2*i+1] = pauli[j]
            creator[2*(i+1)+1] = pauli[j]
            ham.append(''.join(creator))
            coeffs.append(J11)
    return SparsePauliOp(ham, coeffs)

def unitary_time_evolver(ham, *args, num_qbits, time=T, dt=dt):#num_steps=num_time_steps):

    circuit = QuantumCircuit(num_qbits)
    
    for i in range(1, int(time/dt)+1):
        circuit.compose(HamiltonianGate(ham(i*dt, num_qbits, *args), time=dt), inplace=True)
        # print(Operator(HamiltonianGate(ham(i*dt, *args), time=dt)).is_unitary())
    
    return circuit