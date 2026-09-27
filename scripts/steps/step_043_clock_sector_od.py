#!/usr/bin/env python3
"""
Step 043: Clock-Sector Proper-Time Artefact -> Orbit-Determination Rederivation
===============================================================================

PURPOSE:
    Replace the trajectory-perturbation interpretation of the TEP flyby
    response with the canonical clock-sector mechanism of Jakarta v0.14 /
    EFA Yogyakarta: the spacecraft's proper time runs on the matter metric
    g̃, while the orbit-determination (OD) solution assumes the GR time
    standard.  The accumulated proper-time offset along the flyby arc is
    absorbed by the OD least-squares fit and reappears as an *apparent*
    velocity discontinuity — the flyby anomaly.

PIPELINE (Option A rederivation):
    phi(r(t)) -> eta(t) = A(psi) - 1  ->  delta_tau(t) = ∫ eta dt
             -> synthetic clock-carrying injection  z_clock = z_clean + c·eta
                (two-way geometric proxy; not coherent-link clock coupling)
             -> batch least-squares OD fit (GR time standard assumed)
             -> Delta_v_apparent  vs  published anomaly

PHYSICAL MODEL:
    * Excursion field psi(r) = phi_space - phi(r) >= 0 is taken from the
      identical canonical field solution used throughout this pipeline
      (TEPTemporalTopologyModel; surface normalization the Paper 0
      unscreened linear response psi_uns = 1.3926e-9, parameter-free and
      beta-independent; relaxation length lambda_T = 4200 km).  No new
      field model is introduced: the clock channel samples the same
      excursion that the force-channel analyses integrate.
    * Clock-rate offset: eta(t) = exp(beta_A psi/M_Pl) - 1 with the
      universal coupling beta_A = -1  (A < 1 in the well; clocks run
      slower).  This is the AMPLITUDE channel S_A, which is unsuppressed
      at perigee — the shear channel S_Sigma that would screen a force
      does not enter this observable.
    * Observable: the term enters where the link references an onboard
      time standard — a one-way downlink, a non-coherent turnaround, or
      regenerative ranging: a fractional rate offset eta at the
      spacecraft shifts the downlink carrier by eta, mapping to an
      apparent range-rate offset c*eta in the measured Doppler.  The term
      enters ONCE at the transponder — unlike the geometric Doppler it is
      not doubled by the two-way path.  In a fully coherent two-way
      turnaround the transponder multiplies the received carrier phase by
      a fixed ratio with no onboard reference, so the endpoint conformal
      factors cancel identically (nu_ret = R*nu_up) and the term is
      absent — the class in which the published Anderson catalogue was
      recorded (Section 5.2 channel accounting).  The step therefore
      models the clock-carrying observable class: what the excursion
      would reconstruct as on a one-way/non-coherent link over the real
      arcs.  The constant station-side offset
      eta(R_E) is absorbed by DSN clock calibration (Jakarta Rule 8:
      only excursions are observable) and is therefore not modelled.
    * OD filter: identical batch least-squares estimator to Step 012
      (MinimalODFilter, 6-state [x0,v0]; plus the modern empirical-
      acceleration variant for the absorption comparison).  The filter's
      dynamical model is point-mass gravity with NO scalar force — the
      screened shear channel contributes nothing; the only unmodelled
      term is the injected clock signature.

DIFFERENTIAL ESTIMATOR:
    Truth trajectories are the real JPL Horizons state-vector arcs
    (Step 038 output), preserving the actual in/out geometric asymmetry.
    The filter is run twice per flyby — on clean synthetic Doppler and on
    Doppler corrupted by c*eta(t) — and the clock-sector apparent
    velocity shift is the differential
        dv_app = dv_detected(clock-corrupted fit) - dv_detected(clean fit).
    This isolates the clock contribution even where the simplified
    point-mass filter leaves a baseline gravity-field misfit against the
    real arc.

OUTPUT:
    results/step043_clock_sector_od.json — per-flyby eta_perigee,
    accumulated delta_tau, apparent dv (minimal + empirical-accel OD),
    observed anomaly, ensemble correlations, and the fitted global
    clock-response coefficient C_A (the clock-channel analogue of the
    Step 008 beta_fit: a single amplitude normalisation across the
    catalogue, fit by inverse variance).
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.steps.step_007_tep_model import TEPTemporalTopologyModel
from scripts.steps.step_012_od_filter_simulation import (
    MinimalODFilter,
    ModernODFilterWithEmpiricalAccel,
    OrbitalMechanics3D,
    PhysicalConstants,
    SyntheticTrackingNetwork,
)
from scripts.utils.physics import C_LIGHT, M_PL_GEV

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Canonical universal conformal coupling (Jakarta v0.14): A = exp(beta_A phi/M_Pl)
BETA_A = -1.0

# Tracking arc: +/-12 h around perigee captures the approach/departure
# asymmetry; 180 s cadence resolves the lambda_T-scale clock signature.
WINDOW_HOURS = 12.0
CADENCE_S = 180.0

# DSN complexes (lat, lon) in radians — generic three-complex network; the
# per-mission historical tracking schedule is not reconstructed (documented
# limitation; the differential estimator makes the result insensitive to the
# precise network realisation).
DSN_STATIONS = [
    (np.radians(35.4), np.radians(-116.9)),   # Goldstone
    (np.radians(40.4), np.radians(-4.25)),    # Madrid
    (np.radians(-35.4), np.radians(148.98)),  # Canberra
]
TWO_WAY = 2.0
# Filter weighting sigma (0.1 mm/s). Measurements are generated noiselessly
# via compute_range_rate; this sets only the LSQ weight.
NOISE_SIGMA_M_S = 1e-4


class FinePropagator(OrbitalMechanics3D):
    """Point-mass propagator with sub-stepped RK4.

    The base ``propagate`` takes one RK4 step per measurement interval,
    which diverges at 180 s cadence through a ~7000 km perigee pass.
    This subclass subdivides each interval so no RK4 step exceeds
    ``max_step_s``.
    """

    def __init__(self, mu, max_step_s: float = 20.0):
        super().__init__(mu=mu, tep_model=None)
        self.max_step_s = max_step_s

    def propagate(self, initial_state: np.ndarray, t_array: np.ndarray) -> np.ndarray:
        n = len(t_array)
        states = np.zeros((n, 6))
        states[0] = initial_state
        for i in range(1, n):
            dt = t_array[i] - t_array[i - 1]
            n_sub = max(1, int(np.ceil(abs(dt) / self.max_step_s)))
            h = dt / n_sub
            s = states[i - 1]
            for _ in range(n_sub):
                s = self.rk4_step(s, h)
            states[i] = s
        return states


def load_inputs():
    series = json.load(open(RESULTS_DIR / "step038_trajectory_series.json"))["missions"]
    catalog = json.load(open(RESULTS_DIR / "step003_archival_flyby_catalog.json"))["flybys"]
    cat_by_mission = {f["mission_name"]: f for f in catalog}
    return series, cat_by_mission


def extract_arc(mission_series: dict) -> dict:
    """Subsample the real Horizons arc to the +/-WINDOW_HOURS, CADENCE_S window."""
    hours = np.asarray(mission_series["hours_from_perigee"], dtype=float)
    r = np.column_stack([
        np.asarray(mission_series["sc_rx_km"], dtype=float),
        np.asarray(mission_series["sc_ry_km"], dtype=float),
        np.asarray(mission_series["sc_rz_km"], dtype=float),
    ]) * 1e3  # m
    v = np.column_stack([
        np.asarray(mission_series["vx_km_s"], dtype=float),
        np.asarray(mission_series["vy_km_s"], dtype=float),
        np.asarray(mission_series["vz_km_s"], dtype=float),
    ]) * 1e3  # m/s
    mask = np.abs(hours) <= WINDOW_HOURS
    idx = np.where(mask)[0]
    step = max(1, int(round(CADENCE_S / ((hours[1] - hours[0]) * 3600.0))))
    idx = idx[::step]
    t = hours[idx] * 3600.0
    states = np.hstack([r[idx], v[idx]])
    return {"t": t, "states": states}


def clock_signature(model: TEPTemporalTopologyModel, states: np.ndarray) -> np.ndarray:
    """Fractional proper-time rate offset eta = A(psi) - 1 along the arc.

    psi(r) = phi_space - phi(r) >= 0 is the field excursion (v0.14 sign
    convention: phi > 0 near the mass, A < 1 in the well).
    """
    radii = np.linalg.norm(states[:, :3], axis=1)
    psi_gev = np.array([model.phi_space - model.phi(r) for r in radii])
    return np.expm1(BETA_A * psi_gev / M_PL_GEV)


def run_od_fit(t_array, propagator, doppler, measurements, initial_state,
               modern: bool):
    """Run one OD fit; return detected end-minus-start speed change [m/s]."""
    if modern:
        filt = ModernODFilterWithEmpiricalAccel(t_array, propagator, doppler)
        x9 = np.concatenate([initial_state[:6], np.zeros(3)])
        res = filt.run_filter(x9, measurements, t_array)
        states = res["final_states"][:, :6]
    else:
        filt = MinimalODFilter(t_array, propagator, doppler)
        res = filt.estimate(initial_state, measurements)
        states = res["final_states"]
    dv = np.linalg.norm(states[-1, 3:6]) - np.linalg.norm(states[0, 3:6])
    return float(dv), res


def analyse_flyby(name: str, arc: dict, model: TEPTemporalTopologyModel,
                  propagator, network, cat_entry: dict) -> dict:
    t = arc["t"]
    states = arc["states"]
    x0 = states[0].copy()

    eta = clock_signature(model, states)
    delta_tau = np.concatenate([[0.0], np.cumsum(0.5 * (eta[1:] + eta[:-1]) * np.diff(t))])

    # Truth dynamics: point-mass propagation of the real initial state.
    # This is the mechanism-faithful choice, not a numerical convenience:
    # under the clock-sector interpretation the shear channel is screened,
    # so the spacecraft follows GR dynamics and the ONLY non-GR element is
    # the clock rate eta(t), sampled on the real Horizons path (which
    # preserves the physical in/out asymmetry of the excursion).
    truth_states = propagator.propagate(x0, t)
    z_clean = np.array([network.compute_range_rate(truth_states[i], t[i])
                        for i in range(len(t))])
    z_clock = z_clean + C_LIGHT * eta

    dv_min_clean, res_min_clean = run_od_fit(t, propagator, network, z_clean, x0, modern=False)
    dv_min_clock, _ = run_od_fit(t, propagator, network, z_clock, x0, modern=False)
    dv_app_minimal = dv_min_clock - dv_min_clean

    try:
        dv_mod_clean, _ = run_od_fit(t, propagator, network, z_clean, x0, modern=True)
        dv_mod_clock, _ = run_od_fit(t, propagator, network, z_clock, x0, modern=True)
        dv_app_modern = dv_mod_clock - dv_mod_clean
        absorption = 1.0 - dv_app_modern / dv_app_minimal if dv_app_minimal != 0 else np.nan
    except Exception:  # empirical-accel variant is diagnostic-only
        dv_app_modern, absorption = np.nan, np.nan

    # Linear-response cross-check at the clean-fit solution.
    try:
        _, J, _ = MinimalODFilter(t, propagator, network).compute_residuals_and_jacobian(
            res_min_clean["state_estimate"], z_clean)
        dz = C_LIGHT * eta
        dX = np.linalg.lstsq(J.T @ J, J.T @ dz, rcond=None)[0]
        pert = res_min_clean["state_estimate"] + dX
        s_base = propagator.propagate(res_min_clean["state_estimate"], t)
        s_pert = propagator.propagate(pert, t)
        dv_lin = (np.linalg.norm(s_pert[-1, 3:6]) - np.linalg.norm(s_pert[0, 3:6])) - \
                 (np.linalg.norm(s_base[-1, 3:6]) - np.linalg.norm(s_base[0, 3:6]))
    except Exception:
        dv_lin = np.nan

    i_p = int(np.argmin(np.linalg.norm(states[:, :3], axis=1)))
    dtau_in = float(delta_tau[i_p] - delta_tau[0])
    dtau_out = float(delta_tau[-1] - delta_tau[i_p])
    return {
        "mission": name,
        "n_epochs": len(t),
        "r_perigee_km": float(np.linalg.norm(states[i_p, :3]) / 1e3),
        "eta_perigee": float(eta[i_p]),
        "delta_tau_total_s": float(delta_tau[-1] - delta_tau[0]),
        "delta_tau_perigee_s": float(delta_tau[i_p]),
        "delta_tau_inbound_s": dtau_in,
        "delta_tau_outbound_s": dtau_out,
        "dtau_leg_asymmetry_s": dtau_out - dtau_in,
        "dv_app_minimal_od_mm_s": dv_app_minimal * 1e3,
        "dv_app_modern_od_mm_s": (None if np.isnan(dv_app_modern)
                                  else float(dv_app_modern) * 1e3),
        "dv_app_linearized_mm_s": (None if np.isnan(dv_lin) else float(dv_lin) * 1e3),
        "empirical_accel_absorption_fraction": (None if np.isnan(absorption)
                                                else float(absorption)),
        "baseline_postfit_rms_m_s": float(res_min_clean["rms"]),
        "published_anomaly_mm_s": cat_entry.get("published_anomaly_mm_s"),
        "published_anomaly_uncertainty_mm_s": cat_entry.get(
            "published_anomaly_uncertainty_mm_s"),
        "usable_for_analysis": cat_entry.get("usable_for_analysis", False),
    }


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    d = float(np.sqrt((rx @ rx) * (ry @ ry)))
    return float(rx @ ry / d) if d > 0 else np.nan


def main():
    series, cat_by_mission = load_inputs()
    model = TEPTemporalTopologyModel()
    propagator = FinePropagator(mu=PhysicalConstants.MU_EARTH, max_step_s=20.0)
    network = SyntheticTrackingNetwork(
        noise_sigma=NOISE_SIGMA_M_S,
        station_lat_lon_rad=DSN_STATIONS,
        two_way_multiplier=TWO_WAY,
    )

    print(f"Field excursion normalization: phi_space - phi_earth = "
          f"{model.delta_phi:.3e} GeV; eta_surface = "
          f"{np.expm1(BETA_A * model.delta_phi / M_PL_GEV):.3e}")

    rows: list[dict] = []
    for name, mser in series.items():
        if name not in cat_by_mission:
            continue
        cat_entry = cat_by_mission[name]
        if not cat_entry.get("usable_for_analysis", False):
            continue
        arc = extract_arc(mser)
        if len(arc["t"]) < 10:
            continue
        row = analyse_flyby(name, arc, model, propagator, network, cat_entry)
        rows.append(row)
        print(f"  {name:>16}: eta_p={row['eta_perigee']:.2e}  "
              f"dtau={row['delta_tau_total_s']*1e6:.2f} us  "
              f"dv_app={row['dv_app_minimal_od_mm_s']:+.3f} mm/s  "
              f"obs={row['published_anomaly_mm_s']}")

    # Ensemble statistics on the S/N-qualified subset (positive, >3sigma
    # detections are the primary diagnostic rows).
    qual = [r for r in rows if r["published_anomaly_mm_s"] is not None]
    det = [r for r in qual
           if r["published_anomaly_mm_s"] > 0
           and (r["published_anomaly_uncertainty_mm_s"] or 0) > 0
           and r["published_anomaly_mm_s"]
               > 3 * r["published_anomaly_uncertainty_mm_s"]]

    def _corr(subset):
        if len(subset) < 3:
            return {"n": len(subset)}
        app = np.array([r["dv_app_minimal_od_mm_s"] for r in subset])
        obs = np.array([r["published_anomaly_mm_s"] for r in subset])
        return {
            "n": len(subset),
            "pearson_r": float(np.corrcoef(app, obs)[0, 1]),
            "spearman_rho": spearman(app, obs),
        }

    # Global clock-response coefficient C_A: dv_obs = C_A * dv_app
    # (clock-channel analogue of the Step 008 beta_fit amplitude).
    fit_src = det if len(det) >= 3 else qual
    app = np.array([r["dv_app_minimal_od_mm_s"] for r in fit_src])
    obs = np.array([r["published_anomaly_mm_s"] for r in fit_src])
    sig = np.array([max(r["published_anomaly_uncertainty_mm_s"] or 0.05, 0.01)
                    for r in fit_src])
    w = 1.0 / sig**2
    C_A = float((w * obs * app).sum() / (w * app * app).sum()) if np.any(app) else np.nan
    resid = obs - C_A * app
    chi2 = float((w * resid**2).sum())

    output = {
        "step": "step_043_clock_sector_od",
        "description": "Clock-sector proper-time artefact propagated through "
                       "batch least-squares orbit determination (Option A "
                       "rederivation of the flyby response)",
        "model": {
            "field_source": "TEPTemporalTopologyModel canonical excursion "
                            "psi(r) = phi_space - phi(r), lambda_T = 4200 km",
            "beta_A": BETA_A,
            "amplitude_channel": "S_A unsuppressed at perigee; shear channel "
                                 "S_Sigma does not enter the clock observable",
            "observable": "synthetic clock-carrying Doppler proxy; "
                          "z_clock = z_clean + c*eta(t), with a two-way "
                          "geometric baseline and one clock-reference term; "
                          "not a fully coherent two-way link. Constant station "
                          "offset absorbed by calibration",
            "od_filter": "MinimalODFilter (6-state batch LSQ, point-mass "
                         "dynamics, GR time standard) + empirical-acceleration "
                         "variant for absorption diagnostic",
            "estimator": "differential dv_detected(corrupted) - dv_detected(clean)",
            "truth_source": "real JPL Horizons state-vector arcs (step038)",
            "window_hours": WINDOW_HOURS,
            "cadence_s": CADENCE_S,
            "network": "Goldstone/Madrid/Canberra mean, two_way=2, noiseless",
        },
        "field_normalization_gev": {
            "phi_earth": model.phi_earth,
            "phi_space": model.phi_space,
            "delta_phi": model.delta_phi,
            "eta_surface": float(np.expm1(BETA_A * model.delta_phi / M_PL_GEV)),
        },
        "per_flyby": rows,
        "ensemble": {
            "correlations_qualified": _corr(qual),
            "correlations_detections": _corr(det),
            "clock_response_coefficient_C_A": C_A,
            "C_A_fit_subset": [r["mission"] for r in fit_src],
            "C_A_chi2": chi2,
            "C_A_n_dof": len(fit_src) - 1,
            "mean_abs_dv_app_mm_s": float(np.mean(np.abs(app))) if len(app) else None,
        },
    }
    out_path = RESULTS_DIR / "step043_clock_sector_od.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved {out_path}")
    print(f"Ensemble (detections, n={len(det)}): "
          f"{_corr(det)}  C_A={C_A:.3e} chi2={chi2:.1f}")


if __name__ == "__main__":
    main()
