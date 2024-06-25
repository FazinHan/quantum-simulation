Ω = 2.5

h_cut = 1


# In[2]:


import numpy as np

num_layers = 2
num_periods = 1
num_time_steps = 100
T = num_periods * 2*np.pi/Ω
dt = T / num_time_steps


# In[3]:


from qiskit.circuit import ParameterVector, Parameter

cost_threshold = 1e-3


#QFT Circuit
def qft(n):
    """Creates an n-qubit QFT circuit"""
    circuit = QuantumCircuit(n)
    def swap_registers(circuit, n):
        for qubit in range(n//2):
            circuit.swap(qubit, n-qubit-1)
        return circuit
    def qft_rotations(circuit, n):
        """Performs qft on the first n qubits in circuit (without swaps)"""
        if n == 0:
            return circuit
        n -= 1
        circuit.h(n)
        for qubit in range(n):
            circuit.cp(np.pi/2**(n-qubit), qubit, n)
        qft_rotations(circuit, n)
    
    qft_rotations(circuit, n)
    swap_registers(circuit, n)
    return circuit

#Inverse Quantum Fourier Transform
def qft_dagger(qc, n):
    """n-qubit QFTdagger the first n qubits in circ"""
    # Don't forget the Swaps!
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-np.pi/float(2**(j-m)), m, j)
        qc.h(j)
    return qc


def create_ansatz_circuit(qc, num_layers, param_space):
    param_counter = -1
    def ansatz_circuit_0(qc, param_space, param_counter=0):
        print('Number of params:',parameter_space_size)
        # layer 0
        # param_counter=0
        for i in range(qc.num_qubits):
            qc.rx(param_space[param_counter:=param_counter+1],i)
            qc.rz(param_space[param_counter:=param_counter+1],i)
        return param_counter
    def ansatz_circuit_1(qc, param_space, param_counter=0):
        # param_counter = 2 * chain_length
        for i in range(qc.num_qubits-1):
            qc.cx(i,i+1)
        qc.cx(-1,0)
        for i in range(qc.num_qubits):
            qc.rz(param_space[param_counter:=param_counter+1],i)
            qc.rx(param_space[param_counter:=param_counter+1],i)
            qc.rz(param_space[param_counter:=param_counter+1],i)
        return param_counter
    param_counter = ansatz_circuit_0(qc, param_space, param_counter)
    for i in range(num_layers):
        param_counter = ansatz_circuit_1(qc, param_space, param_counter)


def ansatz_circuit_ladder(qc, param_space, layers, entangle_ratio):
    counter = 0
    def layer(qc, params, param_counter):
        for i in range(qc.num_qubits):
            qc.rx(params[param_counter],i)
            param_counter += 1
            qc.rz(params[param_counter],i)
            param_counter += 1
        return param_counter
    def entangle(qc, params, param_counter, double_entangle):
        for i in range(qc.num_qubits//2):
            qc.rzz(params[param_counter], 2*i, 2*i+1)
            param_counter += 1
        if double_entangle:
            for i in range((qc.num_qubits-1)//2):
                qc.rzz(params[param_counter], 2*i+1, 2*i+2)
                param_counter += 1
        return param_counter
    fra = Fraction(entangle_ratio).limit_denominator()
    # print(fra)
    for layer_count in range(layers):
        counter = layer(qc, param_space, counter)
        counter = entangle(qc, param_space, counter, double_entangle=(layer_count%fra.denominator<fra.numerator))
        qc.barrier()



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

def hamiltonian():
    pass


# #### Ladder Hamiltonian

# In[21]:


def entanglement_relator(J, ent_ratio):
    return J*ent_ratio

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

import numpy as np

def calculate_overlaps(ansatz, prev_circuits, parameters, estimator):

    def create_fidelity_circuit(circuit_1, circuit_2):

        """
        Constructs the list of fidelity circuits to be evaluated.
        These circuits represent the state overlap between pairs of input circuits,
        and their construction depends on the fidelity method implementations.
        """
                
        if len(circuit_1.clbits) > 0:
            circuit_1.remove_final_measurements()
        if len(circuit_2.clbits) > 0:
            circuit_2.remove_final_measurements()

        circuit = circuit_1.compose(circuit_2.inverse())
        # circuit.measure_all()
        return circuit
    overlaps = []

    matrix = np.zeros((2**ansatz.num_qubits,2**ansatz.num_qubits))
    matrix[0,0] = 1
    observable = SparsePauliOp.from_operator(matrix)

    for prev_circuit in prev_circuits:
        fidelity_circuit = create_fidelity_circuit(ansatz, prev_circuit)
        estimator_job = estimator.run([(fidelity_circuit, observable, [parameters])])
        estimator_result = estimator_job.result()[0]
        value = estimator_result.data.evs[0]
        
        overlaps.append(value)
    
    return np.array(overlaps)


import numpy as np

def cost_func_vqd(parameters, U_T, ansatz, prev_states, step, betas, estimator, hamiltonian, sign=-1):

    '''
    Estimates <ψ|H|ψ> - λ Σ |<0|(U_θβ†)(U_θ)|0>|²

    Where:
    H = observable
    |ψ> = (U_θ†)(U_T)(U_θ)|0>
    '''

    circuit = ansatz.compose(U_T)
    circuit.compose(ansatz.inverse(),inplace=True)
    estimator_job = estimator.run([(circuit, hamiltonian, [parameters])])

    total_cost = 0

    if step > 1:
        overlaps = calculate_overlaps(ansatz, prev_states, parameters, estimator)
        total_cost = np.sum([np.real(betas[state] * overlap**2) for state, overlap in enumerate(overlaps)])

    estimator_result = estimator_job.result()[0]

    value = estimator_result.data.evs[0] - total_cost

    return value*sign


# ### Linear Driver code (suppressed)

# In[11]:


from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
estimator = Estimator()


def convergence_parameter(ansatz, parameters, U_T):
    circuit = ansatz.assign_parameters(parameters)
    floquet_mode = Statevector.from_instruction(circuit)
    value = floquet_mode.expectation_value(U_T)
    return (1 - np.abs(value))**2


if __name__=="__main__":
    from qiskit.quantum_info import SparsePauliOp, Statevector, Operator, Pauli
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import HamiltonianGate, UGate
    from scipy.optimize import minimize
    from optimparallel import minimize_parallel
    import time
    from fractions import Fraction
    
    # j = 1
    
    chain_length = 4
    
    num_states = 2*chain_length
    entangle_ratio = 1
    
    U_T = unitary_time_evolver(hamiltonian_ladder, entangle_ratio, num_qbits=num_states)
    
    matrix = np.zeros((2**num_states, 2**num_states))
    matrix[0,0] = 1
    observable = SparsePauliOp.from_operator(matrix)
    ground_states = []
    excited_states = [] 
    costs = []
    
    layers = range(5,8)
    
    
    
    t0 = time.perf_counter()
    ti = t0
    for num_layers in layers:
    
         
        param_space2 = ParameterVector('test', 500)
        
        ansatz = QuantumCircuit(num_states)
        # print(num_layers)
        ansatz_circuit_ladder(ansatz, param_space2, num_layers, entangle_ratio)
        parameter_space_size2 = len(ansatz.parameters)
        
        param_space2 = ParameterVector('θ', parameter_space_size2)
        ansatz.assign_parameters(param_space2)
        
        k = 2**num_states
        betas = [5]*k
        x0 = np.zeros(parameter_space_size2)
    
        # print(parameter_space_size2)
        # break
        
        
        
        prev_states = []
        prev_opt_parameters = []
        eigenvalues = []
        ϵ2 = 0
    
    
    # try:
        for step in range(1, k + 1):
            
            result = minimize_parallel(cost_func_vqd, x0, args=(U_T, ansatz, prev_states, step, betas, estimator, observable))#, method="bfgs")
            # print(hamiltonian_ladder(np.pi,4,1))
            
            prev_opt_parameters = result.x
            
    
            ϵ2 += convergence_parameter(ansatz, prev_opt_parameters, U_T)
        
        costs.append(ϵ2**.5)
        ti_new = time.perf_counter()
        print(f'{num_layers}-layer circuit computed in {ti_new-t1}s')
        ti = ti_new
    
    # except Exception as e:
    #     print(e)
    # costs = np.array(costs).reshape(8,2)
    
    t1 = time.perf_counter()
    try:
        name = sys.argv[1]
    except:
        name = 'data'
    
    with open(f'./outputs/{name}.npz','wb') as file:
        np.savez(file, layers=layers, costs=costs)
    
    print('time taken: {:.3f}s'.format(t1-t0))



