def create_ansatz_circuit(qc, num_layers, param_space):
    param_counter = -1
    def ansatz_circuit_0(qc, param_space, param_counter=0):
        for i in range(qc.num_qubits):
            qc.rx(param_space[param_counter:=param_counter+1],i)
            qc.rz(param_space[param_counter:=param_counter+1],i)
        return param_counter
    def ansatz_circuit_1(qc, param_space, param_counter=0):
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


def ansatz_circuit_ladder(qc, param_space, layers, entangle_qubits=[]):
    counter = 0
    def layer(qc, params, param_counter):
        for i in range(qc.num_qubits):
            qc.u(params[param_counter],params[param_counter+1],params[param_counter+2],i)
            param_counter += 3
        return param_counter
    def layer_entangle(qc, params, param_counter):
        for j in range(qc.num_qubits//2):
                qc.rzz(params[param_counter], 2*j, 2*j+1)
                param_counter += 1
        return param_counter
    def entangle(qc, params, param_counter, entangle_qubits_local):
        for i in entangle_qubits_local:
            qc.rzz(params[param_counter], *i)
            param_counter += 1
        return param_counter
    for layer_count in range(layers):
        counter = layer(qc, param_space, counter)
        counter = layer_entangle(qc, param_space, counter)
        # qc.barrier()
    counter = entangle(qc, param_space, counter, entangle_qubits)

if __name__=="__main__":
    from qiskit.circuit import ParameterVector
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(4)
    params = ParameterVector('θ',14*2+2)
    ansatz_circuit_ladder(qc, params, 2, entangle_qubits=[(0,2),(1,3)])
    ansatz_circuit_ladder(qc, params, 2)
    print(qc.draw())