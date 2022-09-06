import numpy as np
from qiskit.circuit import QuantumRegister, QuantumCircuit
from qiskit.circuit.library import MCXGate
from qiskit.extensions import SingleQubitUnitary
from math import log2

def qrd(U):
    dim = len(U)
    n_qubits = (dim).bit_length() - 1
    qubits = QuantumRegister(n_qubits)
    circuit = QuantumCircuit(qubits)
    gates_set = []
    for j in range(dim - 1):
        for i in range(j+1,dim):
            if U[i,j] != 0:
                U_k = compute_U_k(U, i, j)
                U_til = compute_U_til(U_k, i, j)
                gates_set.append(single_qubit_gates([i,j],U_til,n_qubits))
                U = U_k @ U
    for gates in reversed(gates_set):
        cnots_len = len(gates)-1
        for cnot in range(cnots_len):
            circuit.append(gates[cnot]["gate"], gates[cnot]["control_qubits"])
        while gates != []:
            gate = gates.pop()
            circuit.append(gate["gate"], gate["control_qubits"])
    print(circuit)
    return circuit

def compute_U_k(U, i, j):
    dim = len(U)
    U_k = np.eye(dim, dtype="complex_")
    length = np.linalg.norm((U[j,j],U[i,j])) if j != dim - 1 else 1
    U_k[j,j] = np.conj(U[j,j]) / length
    U_k[i,j] = U[i,j] / length
    U_k[j,i] = np.conj(U[i,j]) / length
    U_k[i,i] = -U[j,j] / length
    return U_k

def compute_U_til(U_k, i, j):
    U_til = np.zeros((2,2), dtype="complex_")
    U_til[0,0] = U_k[j,j]
    U_til[1,0] = U_k[i,j]
    U_til[0,1] = U_k[j,i]
    U_til[1,1] = U_k[i,i]
    return U_til.conj().T

def single_qubit_gates(state, U_til, n_qubits):
    controlled_gates = []
    gray_codes = build_gray_codes(state[0], state[1], n_qubits)
    if len(gray_codes) > 2:
        for i in range(len(gray_codes)-2):
            state = [gray_codes[i], gray_codes[i+1]]
            controlled_gates.append(single_qubit_cnot_gate(state, n_qubits))

    state = [gray_codes[-2], gray_codes[-1]]
    controlled_gates.append(single_qubit_controlled_u_gate(state, U_til, n_qubits))
    return controlled_gates

def single_qubit_cnot_gate(state, n_qubits):  
    target_qubit = compute_target_qubit(state[0], state[1])
    control_qubits, control_state = controls(state[0], target_qubit, n_qubits)
    gate = MCXGate(n_qubits-1, ctrl_state=control_state)
    return {'control_qubits': control_qubits, 'gate': gate}

def single_qubit_controlled_u_gate(state, U_til, n_qubits):
    target_qubit = compute_target_qubit(state[0], state[1])
    control_qubits, control_state = controls(state[0], target_qubit, n_qubits)
    gate = SingleQubitUnitary(U_til).control(num_ctrl_qubits=n_qubits-1, ctrl_state=control_state)
    return {'control_qubits': control_qubits, 'gate': gate}

def controls(n, target_qubit, n_qubits):
    control_qubits = []
    control_state = ''
    for i in range(n_qubits):
        bit = 1<<i
        if i != target_qubit:
            control_state = ('1' if n & bit == bit else '0') + control_state
            control_qubits.append(i)
    control_qubits.append(target_qubit)
    return [control_qubits, control_state]

def compute_target_qubit(state1, state2):
    return int(log2(state1 ^ state2))

def hamming_distance(n1, n2):
    return n1 ^ n2

def build_gray_codes(n1, n2, n_qubits):
    if n1 > n2:
        n1, n2 = n2, n1
    codes = [n1]
    for i in range(n_qubits):
        bit = 1<<i
        if (bit & n1 == 0) and (bit & n2 == bit):
            n1 = n1 | bit
        elif (bit & n1 == bit) and (bit & n2 == 0):
            n1 = n1 ^ bit
        if n1 != codes[-1]:
            codes.append(n1)        
    return codes         


