# -*- coding: utf-8 -*-
"""Parameters for the beta-adrenergic signaling model.

This module mirrors the MATLAB ``daePARAMS.m`` file but returns a Python
dictionary instead of a numeric array. Only a subset of parameters frequently
used in the simplified Python translation are provided. Additional values can be
added as needed.
"""
from __future__ import annotations


def dae_params(KR: float, KL: float, KA: float, KG: float, alpha_L: float, alpha_A: float, gamma_L: float, gamma_A: float, factor: float) -> dict[str, float]:
    """Return a dictionary of model parameters."""
    params = {
        "Ltot": 0.0,
        "Atot": 0.0,
        "FSK": 0.0,
        "IBMX": 0.0,
        # Receptor/Gs module
        "b1ARtot": 0.0132,
        "Gstot": 3.83,
        "kf_bARK": 1.1e-6,
        "kr_bARK": 2.2e-6,
        "kf_PKA": 3.6e-6,
        "kr_PKA": 2.2e-6,
        "k_G_act": 16e-3,
        "k_G_hyd": 0.8e-3,
        "k_G_reassoc": 1.21,
        "kr_Ra": KR,
        "kr_LRa1": alpha_L * KL,
        "kr_LRi": KL,
        "kr_RaG": KG,
        "kr_LRaG2": gamma_L * KG,
        "kr_ARa1": alpha_A * KA,
        "kr_ARi": KA,
        "kr_ARaG2": gamma_A * KG,
        # Additional parameters used by ``dae_ode``
        "epsilon": 10.0,
    }
    return params
