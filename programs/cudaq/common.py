import cudaq
from cudaq import spin
import numpy as np
import time
import sys
from typing import List

@cudaq.kernel
def trotter_kernel(n_sites: int, t: float, K: int,
                 coefficients: List[complex], words: List[cudaq.pauli_word], init_angle: float):
    qubits = cudaq.qvector(n_sites)
    for i in range(n_sites):
        ry(init_angle * (1 - 2 * (i % 2)), qubits[i])
    dt = t / K
    for k in range(K):
        for i in range(len(coefficients)):
            exp_pauli(coefficients[i].real * dt, qubits, words[i])

def trotter(n_sites: int, H: cudaq.SpinOperator, t: float, K: int, init_angle: float = 0.0) -> np.ndarray:
    # Extract coefficients and words
    coefficients = [term.evaluate_coefficient() for term in H]
    words = [term.get_pauli_word(H.qubit_count) for term in H]

    # Time evolution
    state = cudaq.get_state(trotter_kernel, n_sites, t, K, coefficients, words, init_angle)
    state = np.conj(np.array(state))
    return state
