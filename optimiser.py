from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qutip import sigmax, sigmay, sigmaz
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from optimparallel import minimize_parallel
import time, sys
import numpy as np
import matplotlib.pyplot as plt
estimator = Estimator()

from physics import hamiltonian_ladder, unitary_time_evolver, T
from ansatzor import ansatz_circuit_ladder
from information import cost_func_vqd, convergence_parameter, determine_next_filename

### Systeme
num_rungs = 1

### Simulatore
num_layers = 1

num_qubits = 2*num_rungs

B = sys.argv[1]

U_T = unitary_time_evolver(hamiltonian_ladder, num_rungs, float(B), num_qbits=num_qubits)

matrix = np.zeros((2**num_qubits, 2**num_qubits))
matrix[0,0] = 1
observable = SparsePauliOp.from_operator(matrix)

singlets = []
triplets = []

layers = [1]

t0 = time.perf_counter()
ti = t0
for num_layers in layers:

     
    param_space2 = ParameterVector('test', 500)
    
    ansatz = QuantumCircuit(num_qubits)
    # print(num_layers)
    ansatz_circuit_ladder(ansatz, param_space2, num_layers)
    parameter_space_size2 = len(ansatz.parameters)
    
    param_space2 = ParameterVector('θ', parameter_space_size2)
    ansatz.assign_parameters(param_space2)
    
    k = 2**num_qubits
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
        

        # ϵ2 += convergence_parameter(ansatz, prev_opt_parameters, U_T)
        floquet_state = Statevector.from_instruction(ansatz.assign_parameters(prev_opt_parameters))
        eigenvalues.append(-np.angle(floquet_state.expectation_value(U_T))/T)
    eigenvalues = np.array(eigenvalues)
    eigenvalues.sort()
    singlets.append(eigenvalues.pop(0))
    for i in eigenvalues:
        triplets.append(i)
    ti_new = time.perf_counter()
    print(f'{num_layers}-layer circuit computed in {ti_new-ti}s')
    ti = ti_new

# except Exception as e:
#     print(e)
# costs = np.array(costs).reshape(8,2)

t1 = time.perf_counter()

folder='data//data1'
with open(determine_next_filename(folder, filetype='npz'),'wb') as file:
    np.savez(file, singlets=singlets, triplets=triplets)

# print('time taken: {:.3f}s'.format(t1-t0))