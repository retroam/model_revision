# -*- coding: utf-8 -*-
"""Python translation of daeODE.m

This module provides a function ``dae_ode`` implementing the beta-adrenergic
signaling model originally written in MATLAB (``daeODE.m``).

The code has been adapted for Python and ``scipy`` integrators. Parameters
should be supplied as a dictionary obtained from :mod:`dae_params`.

The original MATLAB files remain untouched for reference.
"""
from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np


def dae_ode(t: float, y: Iterable[float], params: dict[str, float], flag: int) -> Tuple[np.ndarray, np.ndarray]:
    """Compute derivatives for the signaling model.

    Parameters
    ----------
    t : float
        Current time (ms).
    y : Iterable[float]
        Current state vector.
    params : dict[str, float]
        Model parameters (see :mod:`dae_params`).
    flag : int
        Optional flag used in the original model to control G activation.

    Returns
    -------
    dydt : numpy.ndarray
        Derivatives of the state variables.
    algvars : numpy.ndarray
        Algebraic variables ``[Ra, LRi, LRa, RaG, LRaG, ARi, ARa, ARaG]``.
    """

    y = np.asarray(y)

    # Unpack parameters used below. Only a subset is needed for this simplified
    # translation.
    Ltot = params["Ltot"]
    Atot = params["Atot"]
    b1ARtot = params["b1ARtot"]
    Gstot = params["Gstot"]
    k_G_act = params["k_G_act"]
    k_G_hyd = params["k_G_hyd"]
    k_G_reassoc = params["k_G_reassoc"]
    epsilon = params["epsilon"]

    # Receptor kinetics
    kr_Ra = params["kr_Ra"]
    kr_LRa1 = params["kr_LRa1"]
    kr_LRi = params["kr_LRi"]
    kr_RaG = params["kr_RaG"]
    kr_LRaG2 = params["kr_LRaG2"]
    kr_ARa1 = params["kr_ARa1"]
    kr_ARi = params["kr_ARi"]
    kr_ARaG2 = params["kr_ARaG2"]

    # PKA/BARK rates
    kf_bARK = params["kf_bARK"]
    kr_bARK = params["kr_bARK"]
    kf_PKA = params["kf_PKA"]
    kr_PKA = params["kr_PKA"]

    # Map state vector (see daeODE.m for ordering)
    (
        Ri,
        G,
        b1AR_S464,
        b1AR_S301,
        GsaGTPtot,
        GsaGDP,
        Gsby,
        AC_GsaGTP,
        cAMPtot,
        PDEp,
        RC_I,
        RCcAMP_I,
        RCcAMPcAMP_I,
        RcAMPcAMP_I,
        PKACI,
        PKACI_PKI,
        RC_II,
        RCcAMP_II,
        RCcAMPcAMP_II,
        RcAMPcAMP_II,
        PKACII,
        PKACII_PKI,
        I1p_PP1,
        I1ptot,
        LCCap,
        LCCbp,
        PLBp,
        PLMp,
        TnIp,
    ) = y

    # Extended ternary complex model
    Ra = Ri / kr_Ra
    LRi = Ltot * Ri / kr_LRi
    LRa = Ltot * Ra / kr_LRa1
    RaG = Ra * G / kr_RaG
    LRaG = LRa * G / kr_LRaG2
    ARi = Atot * Ri / kr_ARi
    ARa = Atot * Ra / kr_ARa1
    ARaG = ARa * G / kr_ARaG2
    b1ARact = b1ARtot - b1AR_S464 - b1AR_S301
    dRi = b1ARact - Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG - Ri
    dG = Gstot - LRaG - RaG - ARaG - G

    bARK_desens = kf_bARK * (LRa + LRaG + RaG + ARa + ARaG)
    bARK_resens = kr_bARK * b1AR_S464
    PKA_desens = kf_PKA * PKACI * b1ARact
    PKA_resens = kr_PKA * b1AR_S301
    db1AR_S464 = bARK_desens - bARK_resens
    db1AR_S301 = PKA_desens - PKA_resens

    if flag == 1:
        G_act = k_G_act * (RaG + LRaG + ARaG) * (
            epsilon / (epsilon + Atot)
        )
    else:
        G_act = k_G_act * (RaG + LRaG + ARaG)

    G_hyd = k_G_hyd * GsaGTPtot
    G_reassoc = k_G_reassoc * GsaGDP * Gsby
    dGsaGTPtot = G_act - G_hyd
    dGsaGDP = G_hyd - G_reassoc
    dGsby = G_act - G_reassoc

    # This simplified translation omits the downstream cAMP/PKA equations for
    # brevity. In a full translation these would mirror the MATLAB code.
    dcAMPtot = 0.0
    dAC_GsaGTP = 0.0
    dPDEp = 0.0

    dRC_I = dRCcAMP_I = dRCcAMPcAMP_I = dRcAMPcAMP_I = 0.0
    dPKACI = dPKACI_PKI = 0.0
    dRC_II = dRCcAMP_II = dRCcAMPcAMP_II = dRcAMPcAMP_II = 0.0
    dPKACII = dPKACII_PKI = 0.0
    dI1p_PP1 = dI1ptot = 0.0
    dLCCap = dLCCbp = dPLBp = dPLMp = dTnIp = 0.0

    dydt = np.array(
        [
            dRi,
            dG,
            db1AR_S464,
            db1AR_S301,
            dGsaGTPtot,
            dGsaGDP,
            dGsby,
            dAC_GsaGTP,
            dcAMPtot,
            dPDEp,
            dRC_I,
            dRCcAMP_I,
            dRCcAMPcAMP_I,
            dRcAMPcAMP_I,
            dPKACI,
            dPKACI_PKI,
            dRC_II,
            dRCcAMP_II,
            dRCcAMPcAMP_II,
            dRcAMPcAMP_II,
            dPKACII,
            dPKACII_PKI,
            dI1p_PP1,
            dI1ptot,
            dLCCap,
            dLCCbp,
            dPLBp,
            dPLMp,
            dTnIp,
        ]
    )

    algvars = np.array([Ra, LRi, LRa, RaG, LRaG, ARi, ARa, ARaG])
    return dydt, algvars
