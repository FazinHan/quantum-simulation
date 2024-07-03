from qiskit.quantum_info import SparsePauliOp
import numpy as np
from qiskit.quantum_info import Statevector

def convergence_parameter(ansatz, parameters, U_T):
    circuit = ansatz.assign_parameters(parameters)
    floquet_mode = Statevector.from_instruction(circuit)
    value = floquet_mode.expectation_value(U_T)
    return (1 - np.abs(value))**2

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

def determine_next_filename(filename='output',folder='outputs',filetype='png',exists=False):
    num = 1
    filename = lambda num: f'{filename}{num}.{filetype}'
    import os
    while os.path.isfile(os.path.join('.',folder,filename(num))):
        num += 1
    if exists:
        num -= 1
    return os.path.join('.',folder,filename(num))