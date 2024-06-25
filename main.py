from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.primitives import StatevectorEstimator as Estimator
from qutip import sigmax, sigmay, sigmaz
estimator = Estimator()

chain_length = 1

num_qubits = 2*chain_length