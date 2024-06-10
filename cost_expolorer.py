import numpy as np
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt
from qiskit.circuit import ParameterVector
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qiskit.quantum_info import SparsePauliOp, Statevector, Operator, Pauli
from qiskit.circuit.library import HamiltonianGate, UGate
import os, time

sampler = Sampler()
estimator = Estimator()

chain_length = 4

OO = 2.5

omit_ratio = 0

h_cut = 1

num_layers = 2
num_periods = 1
num_time_steps = 100
T = num_periods * 2*np.pi/OO
dt = T / num_time_steps

def hamiltonian_linear(t, A, Δ=1, omega=OO):
    ham = SparsePauliOp(['Z','X'] , [-Δ/2, A/2*np.cos(omega*t)])
    # plt.plot(t, A*np.cos(OO*t)/2,'.')
    return ham

def unitary_time_evolver(ham, *args, num_qbits=chain_length, time=T, dt=dt):#num_steps=num_time_steps):

    circuit = QuantumCircuit(num_qbits)
    
    for i in range(1, int(time/dt)+1):
        circuit.compose(HamiltonianGate(ham(i*dt, *args), time=dt), inplace=True)
        # print(Operator(HamiltonianGate(ham(i*dt, *args), time=dt)).is_unitary())
    
    return circuit

def calculate_overlaps(ansatz, prev_circuits, parameters, sampler):

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
        circuit.measure_all()
        return circuit
    overlaps = []

    for prev_circuit in prev_circuits:
        fidelity_circuit = create_fidelity_circuit(ansatz, prev_circuit)
        sampler_job = sampler.run([(fidelity_circuit, parameters)])
        meas_data = sampler_job.result()[0].data.meas
        
        counts_0 = meas_data.get_int_counts().get(0, 0)
        shots = meas_data.num_shots
        overlap = counts_0/shots
        overlaps.append(overlap)
    
    return np.array(overlaps)

def cost_func_vqd(parameters, U_T, ansatz, prev_states, step, betas, estimator, sampler, hamiltonian, sign=-1):

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
        overlaps = calculate_overlaps(ansatz, prev_states, parameters, sampler)
        total_cost = np.sum([np.real(betas[state] * overlap**2) for state, overlap in enumerate(overlaps)])

    estimator_result = estimator_job.result()[0]

    value = estimator_result.data.evs[0] - total_cost

    return value*sign

def cost_wrapper(theta, lam):
    ansatz = QuantumCircuit(1)
    thetas = ParameterVector('theta',3)
    ansatz.u(*thetas,0)
    U_T = unitary_time_evolver(hamiltonian_linear, 1/OO, num_qbits=1)
    
    matrix = np.array([[1,0],[0,0]])
    observable = SparsePauliOp.from_operator(matrix)

    return cost_func_vqd(np.array((theta,0,lam)), U_T, ansatz, [], 0, [], estimator, sampler, observable, sign=1)

if __name__=="__main__":
    theta = np.linspace(-2*np.pi,2*np.pi)
    lam = np.linspace(0,2*np.pi)
    theta, lam = np.meshgrid(theta, lam)
    costs = np.zeros(theta.shape)

    t0 = time.perf_counter()
    for i in range(theta.shape[0]):
        for j in range(theta.shape[1]):
            costs[i,j] = cost_wrapper(theta[i,j],lam[i,j])
    t1 = time.perf_counter()
    print('time taken {:.3f} s'.format(t1-t0))

    num = 1
    while os.path.isfile(f'.//outputs//figure{num}.png'):
        num += 1

    plt.contourf(theta, lam, costs)
    plt.savefig(f'.//outputs//figure{num}.png')
    