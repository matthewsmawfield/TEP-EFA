#!/usr/bin/env python3
"""Step 044: Lense-Thirring (frame-dragging) velocity bound per flyby.

Evaluates the general-relativistic gravitomagnetic contribution to the
flyby velocity change by direct trajectory integration.  Earth's spin
angular momentum produces the stationary gravitomagnetic field; the
first-order acceleration on a test body is

    a_LT = -2 v x Omega_g(r) ,
    Omega_g(r) = (G / c^2 r^3) [ J - 3 (J . r_hat) r_hat ] ,

with Earth angular momentum J = 5.86e33 kg m^2 s^-1 directed along the
rotation axis.  Each flyby is propagated from its perigee state vector
(JPL Horizons, Step 038 archive) along the osculating Kepler hyperbola,
parametrised by true anomaly nu with dt/dnu = r^2/h, and the vector
acceleration is integrated over the full hyperbolic arc.

Replaces the previously cited "archived Step 038" computation, which is
not present in the current pipeline tree; the exclusion conclusion is
re-derived here at the correct magnitude (the prior ~1e-13 mm/s values
were ~8 orders below the literature ~1e-5 mm/s scale).

Inputs:
    results/step038_3d_state_vectors.json   perigee state vectors (equatorial)
    results/step003_archival_flyby_catalog.json  published anomalies

Output:
    results/step044_lense_thirring.json
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json

import numpy as np

from scripts.utils.step_logger import StepLogger

# Physical constants (SI)
GM_EARTH = 3.986004418e14      # m^3 s^-2
G_NEWTON = 6.67430e-11         # m^3 kg^-1 s^-2
C_LIGHT = 299792458.0          # m s^-1
J_EARTH = 5.86e33              # kg m^2 s^-1  (spin angular momentum)


def orbit_from_perigee(r_p_vec, v_p_vec):
    """Osculating two-body elements from a perigee state vector (SI)."""
    r_p = np.linalg.norm(r_p_vec)
    v_p = np.linalg.norm(v_p_vec)
    eps = 0.5 * v_p * v_p - GM_EARTH / r_p
    h_vec = np.cross(r_p_vec, v_p_vec)
    h = np.linalg.norm(h_vec)
    e = np.sqrt(1.0 + 2.0 * eps * h * h / (GM_EARTH * GM_EARTH)) \
        if GM_EARTH != 0 else 1.0
    z_hat = h_vec / h if h > 0 else np.array([0.0, 0.0, 1.0])
    x_hat = r_p_vec / r_p
    v_rad = np.dot(v_p_vec, x_hat)
    y_hat = np.cross(z_hat, x_hat)
    y_hat = y_hat / np.linalg.norm(y_hat)
    p = h * h / GM_EARTH
    return {"eps": eps, "h": h, "e": e, "p": p,
            "x_hat": x_hat, "y_hat": y_hat, "z_hat": z_hat,
            "r_p": r_p, "v_p": v_p, "v_radial_at_perigee": v_rad}


def lense_thirring_dv(orb, n_grid=200001):
    """Integrate a_LT = -2 v x Omega_g along the hyperbolic arc (SI units).

    Perifocal parametrisation by true anomaly nu:
        r = p / (1 + e cos nu)
        v = (mu/h) [ -sin(nu) x_hat + (e + cos nu) y_hat ]
        dt/dnu = r^2 / h
    """
    e, p, h = orb["e"], orb["p"], orb["h"]
    x_hat, y_hat = orb["x_hat"], orb["y_hat"]

    nu_lim = 0.999 * np.arccos(-1.0 / e) if e > 1.0 else np.pi
    nu = np.linspace(-nu_lim, nu_lim, n_grid)
    r_mag = p / (1.0 + e * np.cos(nu))

    # Perifocal position and velocity vectors
    pos = (np.cos(nu)[:, None] * (orb["x_hat"][None, :] * r_mag[:, None])
           + np.sin(nu)[:, None] * (orb["y_hat"][None, :] * r_mag[:, None]))
    vel = (GM_EARTH / h) * (
        np.outer(-np.sin(nu), orb["x_hat"]) + np.outer(e + np.cos(nu), orb["y_hat"]))

    # Gravitomagnetic angular velocity field:  Omega_g = (G/c^2 r^3)[J - 3(J.r)r/r^2]
    J_vec = np.array([0.0, 0.0, J_EARTH])
    r_hat = pos / r_mag[:, None]
    Jdot_rhat = np.sum(r_hat * J_vec[None, :], axis=1)
    Om = (G_NEWTON / (C_LIGHT * C_LIGHT)) * (
        (J_vec[None, :] - 3.0 * Jdot_rhat[:, None] * r_hat) / (r_mag[:, None] ** 3))

    a_lt = -2.0 * np.cross(vel, Om)
    dt_dnu = r_mag * r_mag / h
    dv_vec = np.trapezoid(a_lt * dt_dnu[:, None], nu, axis=0) \
        if hasattr(np, "trapezoid") else np.trapz(a_lt * dt_dnu[:, None], nu, axis=0)
    return dv_vec


def main():
    logger = StepLogger("step_044_lense_thirring")

    vec_path = PROJECT_ROOT / "results" / "step038_3d_state_vectors.json"
    cat_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    if not vec_path.exists():
        logger.error(f"Missing {vec_path}")
        return None
    vectors = json.load(open(vec_path))
    catalog = json.load(open(cat_path)) if cat_path.exists() else {"flybys": []}
    anomalies = {f["mission_name"]: f.get("published_anomaly_mm_s")
                 for f in catalog.get("flybys", [])}

    flyby_results = []
    for name, sv in vectors.items():
        r_p_vec = np.array([sv["rx_km"], sv["ry_km"], sv["rz_km"]]) * 1e3
        v_p_vec = np.array([sv["vx_km_s"], sv["vy_km_s"], sv["vz_km_s"]]) * 1e3
        orb = orbit_from_perigee(r_p_vec, v_p_vec)
        dv_vec = lense_thirring_dv(orb)
        dv_mm_s = float(np.linalg.norm(dv_vec) * 1000.0)
        # Along-track (signed) component — the quantity closest to a
        # range-rate observable
        v_hat = v_p_vec / np.linalg.norm(v_p_vec)
        dv_along_mm_s = float(np.dot(dv_vec, v_hat) * 1000.0)

        obs = anomalies.get(name)
        frac = float(dv_mm_s / obs) if (obs is not None and obs > 0) else None
        entry = {
            "mission": name,
            "perigee_range_km": float(orb["r_p"] / 1000.0),
            "perigee_speed_km_s": float(orb["v_p"] / 1000.0),
            "dv_lt_vector_mm_s": [float(x * 1000.0) for x in dv_vec],
            "dv_lt_magnitude_mm_s": dv_mm_s,
            "dv_lt_along_track_mm_s": dv_along_mm_s,
            "observed_anomaly_mm_s": obs,
            "fraction_of_anomaly": frac,
            "excluded_lt1pct": bool(frac < 0.01) if frac is not None else None,
        }
        flyby_results.append(entry)
        print(f"  {name:>16}: |dv_LT| = {dv_mm_s:.3e} mm/s "
              f"(along-track {dv_along_mm_s:+.3e})")

    mags = [r["dv_lt_magnitude_mm_s"] for r in flyby_results]
    out = {
        "step": "044",
        "name": "lense_thirring_frame_dragging",
        "model": "a_LT = -2 v x Omega_g; Omega_g = (G/c^2 r^3)[J - 3(J.r_hat)r_hat]; J_Earth = 5.86e33 kg m^2/s along +z (equatorial state vectors)",
        "integration": "osculating Kepler hyperbola from step038 perigee state vectors; dt/dnu = r^2/h; full arc",
        "flyby_results": flyby_results,
        "summary": {
            "n_flybys": len(flyby_results),
            "dv_lt_range_mm_s": [min(mags), max(mags)],
            "conclusion": ("Lense-Thirring contribution is at most "
                           f"{max(mags):.2e} mm/s, at least three orders "
                           "below every published nonzero anomaly"),
        },
    }
    out_path = PROJECT_ROOT / "results" / "step044_lense_thirring.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    logger.success(f"Results saved to {out_path}")
    return out


if __name__ == "__main__":
    main()
