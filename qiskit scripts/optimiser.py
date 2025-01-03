from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qutip import sigmax, sigmay, sigmaz
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from optimparallel import minimize_parallel
import time
import numpy as np
import matplotlib.pyplot as plt

from physics import hamiltonian_ladder, unitary_time_evolver, T, h_cut
from ansatzor import ansatz_circuit_ladder
from information import cost_func_vqd, convergence_parameter, determine_next_filename, penalty

'''Do not run for more than one layer at a time!'''

estimator = Estimator()

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
        x0 = np.random.uniform(high=1,size=parameter_space_size2)
        # x0 = np.ones(parameter_space_size2)*np.pi/8
    
        # print(parameter_space_size2)
        # break
        
        
        
        prev_states = []
        prev_opt_parameters = []
        eigenvalues = []
        layer_step = []
        costs = []
        ϵ2 = 0
        penalties = []

    
    
    # try:
        for step in range(1, k + 1):
            
            result = minimize_parallel(cost_func_vqd, x0, args=(U_T, ansatz, prev_states, step, betas, estimator, observable))#, method="bfgs")
            
            prev_opt_parameters = result.x
            
            cost = result.fun # np.float64
            
            curr_state = ansatz.assign_parameters(prev_opt_parameters)
            
            prev_states.append(curr_state)
            
            floquet_state = Statevector.from_instruction(curr_state)
            eigenvalues.append([cost, -h_cut*np.angle(floquet_state.expectation_value(U_T))/T]) # append [state_label, floquet_energy]
            
            layer_step.append(step)

            overlap = penalty(prev_opt_parameters, U_T, ansatz, prev_states, step, betas, estimator, observable)
            penalties.append(overlap)
        
        eigenvalues = np.array(eigenvalues)
        penalties = np.array(penalties)
        # singlet = eigenvalues[0]
        # triplets = eigenvalues[1:]
        ti_new = time.perf_counter()
        print(f'{num_layers}-layer circuit computed in {np.round(ti_new-ti, 3)}s')
        ti = ti_new
    return B, eigenvalues, layer_step, penalties