import numpy as np
from qiskit.circuit.library import HamiltonianGate
from qiskit import QuantumCircuit
from qutip import tensor, sigmax, sigmay, sigmaz, qeye

Ω = 2.5
h_cut = 1

num_periods = 1
num_time_steps = 100
T = num_periods * 2*np.pi/Ω
dt = T / num_time_steps

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


def hamiltonian_ladder(t, num_rungs, J=1, ratio=1, B=1, omega=2.5):
    num_states = num_rungs*2
    J11 = J*ratio
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

def qutip_ladder_hamiltonian(num_rungs, J=1, ratio=1, B=1, omega=2.5):
    num_states = num_rungs*2
    J11 = ratio * J
    pauli_list = [qeye(2)]*num_states
    H0 = []
    H1 = []
    def H1_t(t, args):
        return B * np.cos(omega * t)
    H2 = []
    def H2_t(t, args):
        return B * np.sin(omega * t)
    for x in ['x','y','z']:
        for i in range(0,num_states-1,2):
            op = pauli_list[:]
            op[i] = eval(f'sigma{x}()')
            op[i+1] = eval(f'sigma{x}()')
            H0.append(J*tensor(op))
            try:
                op2 = pauli_list[:]
                op2[i] = eval(f'sigma{x}()')
                op2[i+2] = eval(f'sigma{x}()')
                H0.append(J11*tensor(op2))
            except IndexError:
                pass
            try:
                op3 = pauli_list[:]
                op3[i+1] = eval(f'sigma{x}()')
                op3[i+3] = eval(f'sigma{x}()')
                H0.append(J11*tensor(op3))
            except IndexError:
                pass
    for i in range(0,num_states-1,2):   
        op4 = pauli_list[:]
        op4[i] = sigmax()
        op4[i+1] = sigmax()
        op5 = pauli_list[:]
        op5[i] = sigmay()
        op5[i+1] = sigmay()
        H1.append(tensor(op4))
        H2.append(tensor(op5))
    return [*H0, [H1, H1_t], [H2, H2_t]]


def unitary_time_evolver(ham, *args, num_qbits, time=T, dt=dt):#num_steps=num_time_steps):

    circuit = QuantumCircuit(num_qbits)
    
    for i in range(1, int(time/dt)+1):
        circuit.compose(HamiltonianGate(ham(i*dt, num_qbits, *args), time=dt), inplace=True)
        # print(Operator(HamiltonianGate(ham(i*dt, *args), time=dt)).is_unitary())
    
    return circuit

if __name__=="__main__":
    ham = hamiltonian_ladder(0, 4)
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    def update(i):
        ham = hamiltonian_ladder(i/2/np.pi, 4)
        matrice.set_array(ham.to_matrix().real)

    fig, ax = plt.subplots()
    matrice = ax.matshow(ham.to_matrix().real)
    plt.colorbar(matrice)
    
    ani = animation.FuncAnimation(fig, update, frames=150, interval=50*3)
    plt.show()