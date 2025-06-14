# -*- coding: utf-8 -*-
"""Run a simulation of the signaling model using SciPy."""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from dae_params import dae_params
from dae_ode import dae_ode


def run_simulation():
    """Simulate the model for a simple protocol."""
    # Parameter values based on ``daeCompute.m``
    KR = 10.0
    KL = 0.2
    KA = 0.2
    KG = 2.4131
    alpha_L = 1.0
    alpha_A = 1.0
    gamma_L = 0.3762
    gamma_A = 1.0
    scaling_factor = 0.1

    params = dae_params(KR, KL, KA, KG, alpha_L, alpha_A, gamma_L, gamma_A, scaling_factor)

    # Initial conditions (29 state variables, all zeros)
    y0 = np.zeros(29)

    def rhs(t, y):
        dydt, _ = dae_ode(t, y, params, flag=0)
        return dydt

    sol = solve_ivp(rhs, (0, 20 * 60 * 1000), y0, method="LSODA")
    return sol


if __name__ == "__main__":
    sol = run_simulation()
    print("Simulation finished with", len(sol.t), "time points")
