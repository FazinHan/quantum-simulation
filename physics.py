import numpy as np
from qiskit.circuit.library import HamiltonianGate
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qutip import tensor, sigmax, sigmay, sigmaz, qeye, FloquetBasis

'''
All tests show hamiltonians are not equivalent somehow.
'''

num_rungs = 1
num_layers = 1
num_qubits = 2*num_rungs

Ω = 2.5
h_cut = 1

num_periods = 1
num_time_steps = 100
T = num_periods * 2*np.pi/Ω
dt = T / num_time_steps

J = 1
JII = 1
B_range = [0,.2]

def hamiltonian_circular(t, A=2, J=J, omega=Ω):
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

def hamiltonian_linear(t, A=2, Δ=1, omega=Ω):
    ham = SparsePauliOp(['Z','X'] , [-Δ/2, A/2*np.cos(omega*t)])
    # plt.plot(t, A*np.cos(Ω*t)/2,'.')
    qutip_ham = [sigmaz()*-Δ/2,[sigmax(), lambda t, args: args.get('A',A)/2*np.cos(args.get('omega',omega)*t)]]
    return ham, qutip_ham


def hamiltonian_ladder(t, num_rungs, B=2, omega=Ω, J=J, JII=JII):
    '''
    Hamiltonians are created with terms of alternating sides: LRLRLRLR...
    Returns a SparsePauliOp.
    
    For open boundary. Edit the exception handling if you want to include periodic boundary conditions.
    '''
    num_qubits = num_rungs*2
    pauli = ['X','Y','Z']
    creator = ['II']*num_rungs
    ham = []
    coeffs = []
    for r in range(num_rungs):
        for j in pauli:
            plist_I_ = creator[:]
            plist_II_R = creator[:]
            plist_II_L = creator[:]
            
            plist_I_[r] = j+j
            coeffs.append(J/4)
            ham.append(''.join(plist_I_))

            try:
                plist_II_R[r+1] = 'I'+j
                plist_II_R[r] = 'I'+j
                plist_II_L[r+1] = j+'I'
                plist_II_L[r] = j+'I' # open boundary
                coeffs.append(JII)
                coeffs.append(JII)
                ham.append(''.join(plist_II_R))
                ham.append(''.join(plist_II_L))
            except IndexError:
                pass
            # print(plist_I_)
            # print(plist_II_L)
            # print(plist_II_R)
            
        tlist_1 = creator[:]
        tlist_2 = creator[:]
        tlist_3 = creator[:]
        tlist_4 = creator[:]

        tlist_2[r] = 'IX'
        tlist_1[r] = 'XI'
        ham.append(''.join(tlist_1))
        ham.append(''.join(tlist_2))
        coeffs.append(B/2 * np.cos( omega * t))
        coeffs.append(B/2 * np.cos( omega * t))

        tlist_3[r] = 'YI'
        tlist_4[r] = 'IY'
        # print(tlist_1)
        ham.append(''.join(tlist_3))
        ham.append(''.join(tlist_4))
        coeffs.append(B/2 * np.sin( omega * t))
        coeffs.append(B/2 * np.sin( omega * t))

    operator = SparsePauliOp(ham, coeffs)

    # print(operator.to_matrix())
        
    return operator

def qutip_ladder_hamiltonian(num_rungs, B=2, omega=Ω, J=J, JII=JII):
    '''
    Hamiltonians are created with terms of alternating sides: LRLRLRLR...
    Returns in QuTiP readable (time-dependent) hamiltonian list.

    For open boundary. Edit the exception handling if you want to include periodic boundary conditions.
    '''
    num_qubits = num_rungs*2
    pauli_list = [qeye(2)]*num_qubits
    H0 = []
    H1 = []
    def H1_t(t, args): # may have to use 'eval' instead to return with B and omega subsituted?
        return args.get('B',B)/2 * np.cos(args.get('omega',omega) * t)
    H2 = []
    def H2_t(t, args):
        return args.get('B',B)/2 * np.sin(args.get('omega',omega) * t)
    for sigma in [sigmax(), sigmay(), sigmaz()]:
        for r in range(num_rungs):
            op = pauli_list[:]
            op[2*r] = sigma
            op[2*r+1] = sigma
            H0.append(J/4*tensor(op))
            try:
                opL = pauli_list[:]
                opL[2*(r+1)] = sigma
                opL[2*r] = sigma
                H0.append(JII*tensor(opL))
            except IndexError:
                pass
            try:
                opR = pauli_list[:]
                opR[1+2*(r+1)] = sigma
                opR[1+2*r] = sigma # open boundary
                H0.append(JII*tensor(opR))
            except IndexError:
                pass
            # print(op)
            # print(opL)
            # print(opR)
    for i in range(num_qubits):   
        clockwise = pauli_list[:]
        clockwise[i] = sigmax()
        counter = pauli_list[:]
        counter[i] = sigmay()
        # print(clockwise)
        H1.append(tensor(clockwise))
        H2.append(tensor(counter))
    return [*H0, *[[H, H1_t] for H in H1], *[[H, H2_t] for H in H2]]

def unitary_integrator(operator, func, t, args={}):
    if t == 0:
        return qeye(2)
    steps = num_time_steps* ( (t/T)**.5 )
    delta_t = t/steps
    output = 1
    for step in range( int(steps) ):
        term = -1j*func(delta_t*step, args)*operator
        output *= term.expm()
    return output

def unitary(ham_list, T):
    return np.sum([*[(-1j*i*T/h_cut).expm() for i in ham_list if type(i) != list], *[unitary_integrator(i[0], i[1], T) for i in ham_list if type(i) == list]])


def unitary_time_evolver(ham, *args, num_qbits, time=T, dt=dt):#num_steps=num_time_steps):

    circuit = QuantumCircuit(num_qbits)
    
    for i in range(1, int(time/dt)+1):
        circuit.compose(HamiltonianGate(ham(i*dt, *args), time=dt), inplace=True)
        # print(Operator(HamiltonianGate(ham(i*dt, *args), time=dt)).is_unitary())
    
    return circuit

def floquetor(hamiltonian, T, **kwargs): # ditched

    B = kwargs.get('B')
    omega = kwargs.get('omega')
    num_rungs = kwargs.get('rungs')
    
    ham_list = unitary(hamiltonian(num_rungs, B, omega), T)

    # print(ham_list[0])
    
    f_basis = FloquetBasis(ham_list, T)

    return f_basis

### ISOLATED TESTS

if __name__=="__main__1":
    from qiskit.circuit import Parameter
    for j in np.linspace(0,3*np.pi):
        qiskitham = hamiltonian_ladder(j,1)
        ham = qutip_ladder_hamiltonian(1)
        qutipham = np.sum([i if type(i) != list else i[0]*(i[1](j, {})) for i in ham])
        print(qutipham)
        exit()
        if not np.allclose(qutipham.full(), qiskitham.to_matrix()):
            print(qutipham.full())
            print(qiskitham.to_matrix())
        else:
            print(1)

if __name__=="__main__":
   
    import matplotlib.animation as animation
    
    t = 0
    qiskit_ham, qutip_ham_list = hamiltonian_ladder(t,2).to_matrix(), qutip_ladder_hamiltonian(2)
    # qiskit_ham, qutip_ham_list = hamiltonian_linear(t)
    # qiskit_ham = qiskit_ham.to_matrix()
    qutip_ham = np.sum([i if type(i) != list else i[0]*(i[1](t, {})) for i in qutip_ham_list]).full()

    outputs = qiskit_ham.real,qiskit_ham.imag,qutip_ham.real,qutip_ham.imag
    
    print(np.allclose(qutip_ham, qiskit_ham))
    import matplotlib.pyplot as plt
    def update(t):
        # qiskit_ham, qutip_ham_list = hamiltonian_linear(t)
        # qiskit_ham = qiskit_ham.to_matrix()
        qiskit_ham, qutip_ham_list = hamiltonian_ladder(t,2).to_matrix(), qutip_ladder_hamiltonian(2)
        qiskit_ham, qutip_ham = qiskit_ham, np.sum([i if type(i) != list else i[0]*(i[1](t, {})) for i in qutip_ham_list]).full()
        outputs = qiskit_ham.real,qiskit_ham.imag,qutip_ham.real,qutip_ham.imag
        for idx, array in enumerate(matrices):
            for jdx, matrix in enumerate(array):
                matrix.set_array(outputs[2*idx+jdx])
    fig, axs = plt.subplots(2,2)
    matrices = [[axs[i,j].matshow(outputs[2*i+j],cmap='seismic') for j in range(2)] for i in range(2)]
    ani = animation.FuncAnimation(fig, update, frames=150, interval=50*3)
    axs[0,0].set_ylabel('qiskit ham')
    axs[1,0].set_ylabel('qutip ham')
    axs[1,0].set_xlabel('real part')
    axs[1,1].set_xlabel('imaginary part')
    fig.colorbar(matrices[0][0],ax=axs)
    plt.show()

    # print(np.allclose(SparsePauliOp('XY'*6,1).to_matrix(), tensor([tensor(sigmax(),sigmay())]*6).full()))


if __name__=="__main__1":
    qutip_ham_list = qutip_ladder_hamiltonian(1)
    print('\nQUTIP HAMILTONIAN:\n\n------\n')
    qutip_ham = [i.full() for i in qutip_ham_list if type(i)!=list] + [ (i[0]+i[1](0,{})).full() for i in qutip_ham_list if type(i) == list]
    print(np.sum(qutip_ham,axis=0))
    print()

if __name__=="__main__1":
    qiskit_ham = hamiltonian_ladder(0,1)
    print('\nQISKIT HAMILTONIAN:\n\n------\n')
    print(qiskit_ham.to_matrix())
    print()

if __name__=="__main__1": 
    from qiskit.quantum_info import Statevector, Pauli
    time_evo_qiskit = unitary_time_evolver(hamiltonian_ladder,1,num_qbits=2)
    statez = Statevector.from_instruction(time_evo_qiskit)
    from qutip import sesolve, basis
    time_evo_qutip = sesolve(qutip_ladder_hamiltonian(1), tensor([basis(2,0)]*2), np.linspace(0,T,num_time_steps), e_ops=[tensor([sigmax()]*2),tensor([sigmay()]*2),tensor([sigmaz()]*2)])
    if statez.expectation_value(Pauli('X')) == time_evo_qutip.expect[0][-1]:
        for i in range(3):
            plt.plot(tlist, time_evo_qutip.expect[i])
        plt.show()

if __name__=="__main__1":
    unit = unitary(qutip_ladder_hamiltonian(1),np.pi)
    print(unit.eigenenergies().round(5)/T)