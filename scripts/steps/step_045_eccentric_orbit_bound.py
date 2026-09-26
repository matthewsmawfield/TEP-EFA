"""
Step 045: Eccentric-Orbit Clock Bound on Lambda_TEP

Computes the predicted eccentric-orbit clock modulation under the EFA profile
and compares it directly against the Delva et al. (2018) bound from the 
eccentric Galileo satellites GSAT-0201/0202.

Resolves issue 15-11.
"""

import os
import json
import numpy as np
from pathlib import Path
import sys
from scipy.optimize import root_scalar

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scripts.steps.step_007_tep_model import TEPTemporalTopologyModel
from scripts.utils.physics import R_EARTH

BETA_A = -1.0
M_PL_GEV = 1.22e19

# Delva et al. 2018 bounds
DELVA_ALPHA = 0.19e-5
DELVA_ALPHA_ERR = 2.48e-5
DELVA_GR_MODULATION = 5.3e-11
DELVA_MAX_MODULATION = (abs(DELVA_ALPHA) + 2 * DELVA_ALPHA_ERR) * DELVA_GR_MODULATION

# GSAT-0201/0202 orbit parameters (Galileo disposal orbit)
a = 27977.6 * 1000  # semi-major axis in m (approx 28000 km)
e = 0.162  # eccentricity

r_p = a * (1 - e)
r_a = a * (1 + e)

model = TEPTemporalTopologyModel()
original_lambda = model.lambda_tep

def get_eta(r, m):
    psi_gev = m.phi_space - m.phi(r)
    return np.expm1(BETA_A * psi_gev / M_PL_GEV)

def get_delta_eta_for_lambda(lambda_m):
    m = TEPTemporalTopologyModel()
    m.lambda_tep = lambda_m
    eta_p = get_eta(r_p, m)
    eta_a = get_eta(r_a, m)
    return abs(eta_p - eta_a)

eta_p_orig = get_eta(r_p, model)
eta_a_orig = get_eta(r_a, model)
delta_eta_orig = abs(eta_p_orig - eta_a_orig)

# Find maximum lambda_TEP that satisfies the bound
def root_func(lambda_m):
    return get_delta_eta_for_lambda(lambda_m) - DELVA_MAX_MODULATION

# Search for root between 100 km and 4200 km
res = root_scalar(root_func, bracket=[100000, 4200000], method='brentq')
lambda_max = res.root if res.converged else None

results = {
    "delva_bound": {
        "alpha": DELVA_ALPHA,
        "alpha_err": DELVA_ALPHA_ERR,
        "gr_modulation": DELVA_GR_MODULATION,
        "max_anomalous_modulation": DELVA_MAX_MODULATION,
        "reference": "Delva et al. (2018) PRL 121(23), 231101"
    },
    "orbit_gsat_0201": {
        "a_m": a,
        "e": e,
        "r_p_m": r_p,
        "r_a_m": r_a
    },
    "canonical_model": {
        "lambda_tep_m": original_lambda,
        "eta_periapsis": eta_p_orig,
        "eta_apoapsis": eta_a_orig,
        "predicted_modulation": delta_eta_orig,
        "exceeds_bound_by_factor": delta_eta_orig / DELVA_MAX_MODULATION if DELVA_MAX_MODULATION > 0 else float('inf')
    },
    "bounded_model": {
        "lambda_tep_max_m": lambda_max,
        "lambda_tep_max_km": lambda_max / 1000 if lambda_max else None
    },
    "conclusion": f"If the clock-amplitude excursion operates as a lapse-level term, Delva et al. (2018) bounds the exterior vertical relaxation scale to lambda_TEP <= {lambda_max / 1000:.1f} km, decoupling it from the lateral correlation scale L_c = 4201 km. Alternatively, the excursion enters observables only through non-integrable transport (the holonomy sector), meaning Step 043 uses the wrong observable for this channel."
}

# Ensure results directory exists
out_dir = Path(__file__).resolve().parent.parent.parent / "results"
out_dir.mkdir(exist_ok=True)
out_file = out_dir / "step_045_eccentric_orbit_bound.json"

with open(out_file, "w") as f:
    json.dump(results, f, indent=4)

print(f"Results saved to {out_file}")
print(json.dumps(results, indent=2))
