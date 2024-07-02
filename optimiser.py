from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qutip import sigmax, sigmay, sigmaz
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from optimparallel import minimize_parallel
import time
import numpy as np
import matplotlib.pyplot as plt
estimator = Estimator()

from physics import hamiltonian_ladder, unitary_time_evolver, T, h_cut
from ansatzor import ansatz_circuit_ladder
from information import cost_func_vqd, convergence_parameter, determine_next_filename

def optimiser_main(B, num_rungs = 1, layers = [1]):
    ### Systeme
    
    
    ### Simulatore
    
    
    num_qubits = 2*num_rungs
    
    U_T = unitary_time_evolver(hamiltonian_ladder, num_rungs, B, num_qbits=num_qubits)
    
    matrix = np.zeros((2**num_qubits, 2**num_qubits))
    matrix[0,0] = 1
    observable = SparsePauliOp.from_operator(matrix)
    
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
        costs = []
        layer_step = []
        ϵ2 = 0
    
    
    # try:
        for step in range(1, k + 1):
            
            result = minimize_parallel(cost_func_vqd, x0, args=(U_T, ansatz, prev_states, step, betas, estimator, observable))#, method="bfgs")
            # print(hamiltonian_ladder(np.pi,4,1))
            
            prev_opt_parameters = result.x
            
            cost = result.fun
            
    
            # ϵ2 += convergence_parameter(ansatz, prev_opt_parameters, U_T)
            floquet_state = Statevector.from_instruction(ansatz.assign_parameters(prev_opt_parameters))
            eigenvalues.append(-h_cut*np.angle(floquet_state.expectation_value(U_T))/T)
            costs.append(cost)
            layer_step.append([num_layers, steps])
        eigenvalues = np.array(eigenvalues)
        eigenvalues.sort()
        singlet = eigenvalues[0]
        triplets = eigenvalues[1:]
        ti_new = time.perf_counter()
        costs = np.array(costs)
        layer_step = np.array(layer_step)
        print(f'{num_layers}-layer circuit computed in {ti_new-ti}s')
        ti = ti_new
    return B, singlet, triplets, costs, layer_step