"""Deterministic geometry envelope for TEP flyby predictions (path-resolved; Step 039 uses Step 008 pooled scale)."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from scripts.utils.physics import (
    C_LIGHT,
    DISFORMAL_COUPLING_STRENGTH,
    DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    J2_EARTH,
    J3_EARTH,
    J4_EARTH,
    KG_M3_TO_GEV4,
    LAMBDA_TEP_M,
    M_EARTH,
    M_PL_GEV,
    R_EARTH,
)

# |Ξ| bands for regime labels (manuscript-aligned; thresholds are documented outputs).
XI_REGIME_LOW: float = 0.03
XI_REGIME_HIGH: float = 0.25

# Near-symmetric trajectories require stronger multipole cancellation than the
# leading-order cos(delta_in) - cos(delta_out) factor alone.
ASYMMETRY_COHERENCE_THRESHOLD = 0.12

# Nominal geometry-envelope heuristics (Step 007 / Step 041 stress tests).
# Each carries a nominal ±50% systematic band in Monte Carlo sensitivity sweeps;
# they are not least-squares free parameters in Step 008 (k_fit = 1 for TEP restricted).
GEOMETRY_ENVELOPE_HEURISTIC_KEYS = (
    "inclination_scale",
    "j2_cos2_amp",
    "j2_alt_decay_km",
    "plasma_density_scale_cm3",
    "plasma_density_exponent",
    "velocity_screening_exponent",
    "asymmetry_coherence_threshold",
)


def default_geometry_envelope_heuristics() -> dict[str, float]:
    """Return a fresh copy of nominal envelope heuristics (Step 007 defaults)."""
    return {
        "inclination_scale": 0.15,
        "j2_cos2_amp": 0.00054,
        "j2_alt_decay_km": 2000.0,
        "plasma_density_scale_cm3": 5000.0,
        "plasma_density_exponent": -0.3,
        "velocity_screening_exponent": 4.0,
        "asymmetry_coherence_threshold": float(ASYMMETRY_COHERENCE_THRESHOLD),
    }


def geometry_modulation_factors_from_heuristics(
    altitude_km: float,
    latitude_deg: float,
    velocity_km_s: float,
    plasma_density_cm3: float,
    heuristics: dict[str, float],
    *,
    v_threshold_km_s: float = DISFORMAL_VELOCITY_THRESHOLD_KM_S,
) -> dict[str, float]:
    """
    Deterministic modulation factors for the geometry envelope core (no plasma S_ansatz).

    ``heuristics`` must contain all keys in ``GEOMETRY_ENVELOPE_HEURISTIC_KEYS``.
    """
    missing = [k for k in GEOMETRY_ENVELOPE_HEURISTIC_KEYS if k not in heuristics]
    if missing:
        raise ValueError(
            f"geometry envelope heuristics missing keys: {missing}; "
            f"required: {GEOMETRY_ENVELOPE_HEURISTIC_KEYS}"
        )
    h = heuristics
    f_inclination = 1.0 + h["inclination_scale"] * abs(np.sin(np.radians(latitude_deg)))
    f_j2 = (1.0 - h["j2_cos2_amp"] * np.cos(np.radians(latitude_deg)) ** 2) * np.exp(
        -altitude_km / h["j2_alt_decay_km"]
    )
    f_plasma = (
        1.0 + plasma_density_cm3 / h["plasma_density_scale_cm3"]
    ) ** h["plasma_density_exponent"]
    exp_v = h["velocity_screening_exponent"]
    f_velocity = (
        (v_threshold_km_s / velocity_km_s) ** exp_v
        if velocity_km_s > v_threshold_km_s
        else 1.0
    )
    f_total_core = f_inclination * f_j2 * f_plasma * f_velocity
    return {
        "f_inclination": float(f_inclination),
        "f_j2": float(f_j2),
        "f_plasma": float(f_plasma),
        "f_velocity": float(f_velocity),
        "f_total_core": float(f_total_core),
    }


def zonal_harmonic_bracket(latitude_deg: float, r_m: float) -> float:
    """J2/J3/J4 zonal bracket scaled by (R_Earth / r)^2."""
    sin_lat = math.sin(math.radians(latitude_deg))
    r_ratio_sq = (R_EARTH / r_m) ** 2
    p4 = (35.0 * sin_lat**4 - 30.0 * sin_lat**2 + 3.0) / 8.0
    return (J2_EARTH + J3_EARTH * sin_lat + J4_EARTH * p4) * r_ratio_sq


def j2_only_bracket(latitude_deg: float, r_m: float) -> float:
    """Legacy J2-only bracket retained for component decomposition."""
    sin_lat = math.sin(math.radians(latitude_deg))
    r_ratio_sq = (R_EARTH / r_m) ** 2
    return J2_EARTH * r_ratio_sq


def near_symmetry_cancellation_factor(
    cos_asymmetry: float,
    coherence_threshold: float = ASYMMETRY_COHERENCE_THRESHOLD,
) -> float:
    """
    Quadratic suppression for trajectories with |cos asymmetry| below the
    coherence threshold. Unity for asymmetric detections.
    """
    if coherence_threshold <= 0:
        raise ValueError("coherence_threshold must be positive")
    asym = abs(cos_asymmetry)
    if asym >= coherence_threshold:
        return 1.0
    return (asym / coherence_threshold) ** 3


def disformal_envelope_factor(
    v_sc_m_s: float,
    cos_asymmetry: float,
    v_cmb_frame_kms: float | None = None,
    sc_cmb_cos_theta: float | None = None,
) -> float:
    """
  Disformal response with asymmetry-gated amplitude. The conformal-gradient
  term already carries cos asymmetry; only the disformal excess is gated.
  """
    v_th = DISFORMAL_VELOCITY_THRESHOLD_KM_S
    alpha_b = DISFORMAL_COUPLING_STRENGTH
    asym = abs(cos_asymmetry)

    if v_cmb_frame_kms is not None and sc_cmb_cos_theta is not None and sc_cmb_cos_theta > 0:
        v_km_s = v_cmb_frame_kms
    else:
        v_km_s = v_sc_m_s / 1e3

    if v_km_s > v_th and cos_asymmetry < 0:
        raw = -1.0 * abs(1.0 + alpha_b * (v_km_s / v_th))
    else:
        raw = 1.0 + alpha_b * (v_km_s / v_th) * math.copysign(1.0, cos_asymmetry or 0.0)

    return 1.0 + (raw - 1.0) * asym


def compute_disformal_transition_xi(
    *,
    velocity_km_s: float,
    cos_asymmetry: float,
    field_gradient_ratio: float,
    v_trans_km_s: float = DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    xi_low: float = XI_REGIME_LOW,
    xi_high: float = XI_REGIME_HIGH,
) -> dict[str, Any]:
    """
    Disformal transition scalar Ξ (velocity-activated, manuscript definition).

    Ξ = (v / v_trans)² × |asym| × (|∇φ| / |∇φ_⊕|) × sgn(asym)

    ``field_gradient_ratio`` must be the ratio of field-gradient magnitudes at
    perigee to the reference surface value (|∇φ(r_p)| / |∇φ(R_⊕)|), both from the
    same screening model instance.
    """
    if v_trans_km_s <= 0:
        raise ValueError("v_trans_km_s must be positive")
    if field_gradient_ratio < 0:
        raise ValueError("field_gradient_ratio must be non-negative")

    v_ratio = float(velocity_km_s) / float(v_trans_km_s)
    asym = float(cos_asymmetry)
    sgn = math.copysign(1.0, asym) if asym != 0.0 else 0.0
    xi = (v_ratio**2) * abs(asym) * float(field_gradient_ratio) * sgn
    xi_abs = abs(xi)

    if xi_abs < xi_low:
        regime = "conformal_dominated"
    elif xi_abs < xi_high:
        regime = "mixed"
    else:
        regime = "disformal_dominated"

    return {
        "xi": float(xi),
        "xi_abs": float(xi_abs),
        "regime_classification": regime,
        "v_trans_km_s": float(v_trans_km_s),
        "v_over_v_trans": float(v_ratio),
        "cos_dec_asymmetry": float(asym),
        "field_gradient_ratio": float(field_gradient_ratio),
        "definition": "(v/v_trans)^2 * |asym| * (|grad_phi|/|grad_phi_earth|) * sgn(asym)",
    }


def derive_disformal_transition_scale(
    *,
    lambda_m: float = LAMBDA_TEP_M,
    radius_m: float = R_EARTH,
    mean_density_kg_m3: float | None = None,
    beta_A: float = -1.0,
    epsilon_phi_calibrated: float = 1.0e-8,
    c_m_s: float = C_LIGHT,
) -> dict[str, Any]:
    """
    Audit of the disformal transition-velocity claim (manuscript §3.5).

    The balance condition with the velocity contraction made explicit is

        B(phi) (v . grad phi)^2 / c^2 ~ A^2(phi) - 1.

    Canonical-action assessment: on the solved uniform-sphere Helmholtz
    profile, the microscopic envelope tail B(phi) ~ B_0 varphi^2 at the
    Earth-surface field would require B_0 ~ 1.6e18 to produce a km/s
    transition -- ~21 orders above the calibrated Paper-0/28 envelope
    normalization and ~26 orders above the projected 1e-19-fractional
    holonomy bound (|B_0| <= 5e-8). No kilometre-per-second disformal
    transition therefore exists under the canonical action, and v_trans is
    carried by the pipeline as an empirical response-template scale
    (phenomenological ansatz, ±20%), not as a field-equation output.

    The earlier coefficient B_eff = 4|beta_A|R/lambda is likewise disclosed
    as a bookkeeping identity: defined through the balance itself via
    B_eff = (A^2 - 1)/(v_trans/c)^2, its insertion returns the assumed scale
    (b_eff_implied == b_eff_shell_aspect identically). The closed form

        v_trans = (c/sqrt(2)) * sqrt(lambda/R) * sqrt(epsilon_phi),

    with epsilon_phi = |grad phi_surface| lambda / M_Pl, is retained as a
    dimensional anchor only: evaluated at the UCD-pinned surface combination
    it numerically brackets the template value without deriving it.

    Two evaluations of epsilon_phi are reported:

    - ``epsilon_phi_solved``: linear uniform-sphere Helmholtz response at the
      mean Earth density (underestimates the calibrated amplitude because the
      UCD-saturated interior pins the field above the linear response);
    - ``epsilon_phi_calibrated``: UCD/GNSS-calibrated surface field
      combination (~1e-8).

    Returns a JSON-serialisable dict so the template scale used by the
    pipeline is traceable to evaluated inputs and an explicit claim status.
    """
    if lambda_m <= 0 or radius_m <= 0:
        raise ValueError("lambda_m and radius_m must be positive")
    if mean_density_kg_m3 is None:
        mean_density_kg_m3 = M_EARTH / (4.0 / 3.0 * math.pi * R_EARTH**3)

    hbarc_gev_m = 1.97327e-16  # GeV m (reduced Planck constant x c)
    lam_gevinv = lambda_m / hbarc_gev_m
    rho_gev4 = mean_density_kg_m3 * KG_M3_TO_GEV4

    x = radius_m / lambda_m
    sinh_x = math.sinh(x)
    cosh_x = math.cosh(x)
    e_x = math.exp(x)

    # Uniform-sphere Helmholtz solution: lap phi - phi/lam^2 = -(beta_A/M_Pl) rho
    # Particular solution phi_p = beta_A rho lam^2 / M_Pl (dimensionless varphi_p).
    varphi_p = abs(beta_A) * rho_gev4 * lam_gevinv**2 / M_PL_GEV**2

    # Matched interior/exterior (sinh(r/lam)/r and e^{-r/lam}/r) at r = R.
    varphi_surface = varphi_p * (1.0 - (x + 1.0) * sinh_x / (x * e_x))
    grad_varphi_surface = (
        varphi_p * (x + 1.0) * (x * cosh_x - sinh_x) / (x * e_x * radius_m)
    )
    epsilon_phi_solved = grad_varphi_surface * lambda_m

    geom = math.sqrt(lambda_m / radius_m)
    v_solved_m_s = c_m_s / math.sqrt(2.0) * geom * math.sqrt(epsilon_phi_solved)
    v_calibrated_m_s = (
        c_m_s / math.sqrt(2.0) * geom * math.sqrt(epsilon_phi_calibrated)
    )

    # UCD saturation-pinned evaluation: inside R_sol the field is pinned at the
    # saturation-density source, so the surface field combination is the
    # particular solution evaluated at rho_T rather than the mean-density
    # linear response. epsilon_phi_ucd = |beta_A| rho_T lam^2 / M_Pl^2.
    RHO_T_KG_M3 = 20000.0  # 20 g/cm^3 (GNSS-calibrated, physics.py RHO_T; +-40%)
    rho_T_gev4 = RHO_T_KG_M3 * KG_M3_TO_GEV4
    epsilon_phi_ucd = abs(beta_A) * rho_T_gev4 * lam_gevinv**2 / M_PL_GEV**2
    v_ucd_m_s = c_m_s / math.sqrt(2.0) * geom * math.sqrt(epsilon_phi_ucd)
    # +40% rho_T upper edge of the documented GNSS uncertainty band
    epsilon_phi_ucd_hi = epsilon_phi_ucd * 1.4
    v_ucd_hi_m_s = c_m_s / math.sqrt(2.0) * geom * math.sqrt(epsilon_phi_ucd_hi)

    # Implied shell coefficient: B_eff = (A^2 - 1)/(v/c)^2 at the calibrated
    # scale, with A^2 - 1 ~ 2|beta_A| epsilon_phi. Independent of epsilon_phi:
    # B_eff = 4|beta_A| R / lambda.
    a2_minus_1_calibrated = 2.0 * abs(beta_A) * epsilon_phi_calibrated
    b_eff_implied = a2_minus_1_calibrated / (v_calibrated_m_s / c_m_s) ** 2
    b_eff_shell_aspect = 4.0 * abs(beta_A) * radius_m / lambda_m

    # Literal reading diagnostic: microscopic envelope B(varphi) = B_0 varphi^2
    # at the solved surface field. The B_0 the literal balance would require at
    # v = v_trans shows that the microscopic tail cannot set the scale.
    a2_minus_1_surface = 2.0 * abs(beta_A) * varphi_surface
    b0_required_literal = a2_minus_1_surface / (
        varphi_surface**2 * (v_calibrated_m_s / c_m_s) ** 2
    )

    return {
        "claim_status": "phenomenological_template_scale",
        "status_note": "v_trans is an empirical response-template scale (±20%), not a field-equation output; under the canonical B(phi) envelope no km/s disformal transition exists (see literal_microscopic_balance). The closed form below is retained as a dimensional anchor only.",
        "formula": "v_trans = (c/sqrt(2)) * sqrt(lambda/R) * sqrt(epsilon_phi)",
        "balance_condition": "B_eff (v/c)^2 = A^2 - 1; B_eff = (A^2 - 1)/(v_trans/c)^2 is a bookkeeping identity of the template (b_eff_implied == b_eff_shell_aspect identically), not an action-level coefficient",
        "inputs": {
            "lambda_m": float(lambda_m),
            "radius_m": float(radius_m),
            "mean_density_kg_m3": float(mean_density_kg_m3),
            "beta_A": float(beta_A),
            "x_radius_over_lambda": float(x),
        },
        "helmholtz_solved": {
            "varphi_particular": float(varphi_p),
            "varphi_surface": float(varphi_surface),
            "epsilon_phi_solved": float(epsilon_phi_solved),
            "v_trans_solved_km_s": float(v_solved_m_s / 1e3),
            "note": "Linear uniform-sphere response; underestimates the UCD-calibrated surface field combination.",
        },
        "ucd_saturation_pinned": {
            "rho_T_kg_m3": RHO_T_KG_M3,
            "epsilon_phi_ucd": float(epsilon_phi_ucd),
            "derivation": "epsilon_phi = |beta_A| rho_T lambda^2 / M_Pl^2 -- field pinned at the UCD saturation density inside R_sol, not the mean-density linear response",
            "v_trans_ucd_km_s": float(v_ucd_m_s / 1e3),
            "v_trans_ucd_plus_40pct_rhoT_km_s": float(v_ucd_hi_m_s / 1e3),
            "note": "Central rho_T gives ~14 km/s; the +40% GNSS uncertainty edge gives ~16.5 km/s, bracketing the 16.8 km/s value used by the pipeline.",
        },
        "calibrated": {
            "epsilon_phi_calibrated": float(epsilon_phi_calibrated),
            "epsilon_phi_source": "UCD saturation-pinned surface field combination (see ucd_saturation_pinned)",
            "v_trans_calibrated_km_s": float(v_calibrated_m_s / 1e3),
        },
        "v_trans_used_km_s": float(DISFORMAL_VELOCITY_THRESHOLD_KM_S),
        "b_eff_implied": float(b_eff_implied),
        "b_eff_shell_aspect": float(b_eff_shell_aspect),
        "literal_microscopic_balance": {
            "b0_required_for_v_trans": float(b0_required_literal),
            "note": "B_0 required if the canonical weak-field envelope tail at the solved surface field set the scale; exceeds the calibrated Paper-0/28 envelope normalization by ~21 orders and the projected 1e-19-fractional holonomy bound (|B_0| <= 5e-8) by ~26 orders. No km/s disformal transition therefore exists under the canonical action; the fitted b_disf consistent with zero is the expected admissible-branch outcome.",
        },
    }


def compose_geometry_envelope(
    *,
    altitude_km: float,
    latitude_deg: float,
    velocity_km_s: float,
    cos_asymmetry: float,
    plasma_density_cm3: float,
    modulation: dict[str, float],
    plasma_screening_factor: float,
    plasma_sign_factor: float,
    asymmetry_coherence_threshold: float = ASYMMETRY_COHERENCE_THRESHOLD,
) -> dict[str, Any]:
    """Aggregate deterministic geometry, plasma, and symmetry factors."""
    asymmetry_factor = near_symmetry_cancellation_factor(
        cos_asymmetry, coherence_threshold=asymmetry_coherence_threshold
    )
    envelope = (
        modulation["f_inclination"]
        * modulation["f_j2"]
        * modulation["f_plasma"]
        * modulation["f_velocity"]
        * asymmetry_factor
    )
    return {
        "f_inclination": modulation["f_inclination"],
        "f_j2": modulation["f_j2"],
        "f_plasma_core": modulation["f_plasma"],
        "f_velocity": modulation["f_velocity"],
        "f_geometry_core": modulation["f_total_core"],
        "plasma_density_cm3": plasma_density_cm3,
        "plasma_screening_factor": plasma_screening_factor,
        "plasma_sign_factor": plasma_sign_factor,
        "asymmetry_cancellation_factor": asymmetry_factor,
        "geometry_envelope": envelope,
    }
