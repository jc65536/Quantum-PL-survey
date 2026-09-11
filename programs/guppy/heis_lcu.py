from __future__ import annotations

from typing import Any, Dict

import numpy as np

from . import lcu_common
from ..common import pauli_models


def run_simulation(config: Dict[str, Any]) -> np.ndarray:
    num_sites = int(config["num_sites"])
    total_time = float(config["time"])
    no_loop = bool(config.get("no_loop", True))
    params = config.get("params", {})
    J = float(params.get("J", 1.0))
    field = float(params.get("field", 0.2))

    ham = pauli_models.heisenberg_pauli_terms(num_sites, J, field)
    gamma = pauli_models.taylor_coefficients(ham, total_time)
    coeffs, paulis, phases = pauli_models.lcu_weights_from_gamma(gamma)

    state = lcu_common.lcu_state(coeffs, paulis, phases, no_loop, init_angle=0.25)
    return state


if __name__ == "__main__":
    cfg = {
        "num_sites": 3,
        "time": 0.3,
        "params": {"J": 1.0, "field": 0.7, "trotter_steps": 48},
    }
    state = run_simulation(cfg)
    print(state)
