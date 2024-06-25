from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qutip import sigmax, sigmay, sigmaz
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from optimparallel import minimize_parallel
import time
import numpy as np
estimator = Estimator()

from physics import hamiltonian_ladder, unitary_time_evolver
from ansatzor import ansatz_circuit_ladder
from information import cost_func_vqd, convergence_parameter

### Systeme
chain_length = 1

### Simulatore
num_layers = 1

num_qubits = 2*chain_length

U_T = unitary_time_evolver(hamiltonian_ladder, num_qbits=num_qubits)

matrix = np.zeros((2**num_qubits, 2**num_qubits))
matrix[0,0] = 1
observable = SparsePauliOp.from_operator(matrix)
ground_states = []
excited_states = [] 
costs = []

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
        

        ϵ2 += convergence_parameter(ansatz, prev_opt_parameters, U_T)
    
    costs.append(ϵ2**.5)
    ti_new = time.perf_counter()
    print(f'{num_layers}-layer circuit computed in {ti_new-ti}s')
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
    np.savez(file, layers=layers, costs=costs, params=prev_opt_parameters)

print('time taken: {:.3f}s'.format(t1-t0))