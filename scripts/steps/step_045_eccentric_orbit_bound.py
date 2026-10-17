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
from scripts.utils.physics import R_EARTH, M_PL_GEV
from core import constants as tep_const
from core.scalar_field import screening_y_two_branch

BETA_A = -1.0

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

# Search for root between 100 km and 4200 km (unscreened Yukawa proxy diagnostic)
res = root_scalar(root_func, bracket=[100000, 4200000], method='brentq')
lambda_max = res.root if res.converged else None

# Master-action two-branch kinetic screening evaluation (Paper 0 v0.15):
# du/dr = y_tb(g(r)) * g(r) / c^2 across [r_p, r_a]
r_grid = np.linspace(r_p, r_a, 2000)
g_grid = tep_const.G_NEWTON * tep_const.M_EARTH / (r_grid ** 2)
y_tb_grid = screening_y_two_branch(g_grid)
du_master = float(np.trapezoid(y_tb_grid * g_grid / (tep_const.C_LIGHT ** 2), r_grid))
dphi_gr = float(tep_const.G_NEWTON * tep_const.M_EARTH / (tep_const.C_LIGHT ** 2) * (1.0 / r_p - 1.0 / r_a))
alpha_master = du_master / dphi_gr
delva_alpha_max_2sig = abs(DELVA_ALPHA) + 2.0 * DELVA_ALPHA_ERR
safety_margin_master = delva_alpha_max_2sig / alpha_master

results = {
    "delva_bound": {
        "alpha": DELVA_ALPHA,
        "alpha_err": DELVA_ALPHA_ERR,
        "alpha_2sigma_upper": delva_alpha_max_2sig,
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
    "master_action_two_branch": {
        "kinetic_sector": "selected two-branch P_,X = k sqrt(xi) + 2 xi (k = 16.03, zero free parameters)",
        "dPhi_GR_over_c2": dphi_gr,
        "predicted_modulation_du": du_master,
        "predicted_alpha": alpha_master,
        "delva_2sigma_bound_alpha": delva_alpha_max_2sig,
        "passes_delva_bound": bool(alpha_master < delva_alpha_max_2sig),
        "safety_margin_factor": safety_margin_master
    },
    "canonical_model": {
        "lambda_tep_m": original_lambda,
        "eta_periapsis": eta_p_orig,
        "eta_apoapsis": eta_a_orig,
        "predicted_modulation": delta_eta_orig,
        "exceeds_bound_by_factor": delta_eta_orig / DELVA_MAX_MODULATION if DELVA_MAX_MODULATION > 0 else float('inf'),
        "note": "Unscreened Yukawa proxy (y=1); omits the master-action kinetic screening factor y_tb(g) ~ 7.5e-7 at r ~ 28,000 km."
    },
    "bounded_model": {
        "lambda_tep_max_m": lambda_max,
        "lambda_tep_max_km": lambda_max / 1000 if lambda_max else None
    },
    "conclusion": (
        f"Under the Paper 0 selected two-branch master action (du/dr = y_tb(g) g/c^2), "
        f"kinetic screening at r = 23,445-32,510 km (y_tb ~ 7.5e-7) suppresses the "
        f"GSAT-0201/0202 eccentric-orbit clock modulation to delta_eta = {du_master:.2e} "
        f"(alpha = {alpha_master:.2e}), passing the Delva et al. (2018) 2-sigma bound "
        f"(alpha < {delva_alpha_max_2sig:.2e}) by a factor of {safety_margin_master:.1f}x "
        f"with zero free parameters (whereas an unscreened y=1 Yukawa profile would "
        f"require lambda_TEP <= {lambda_max / 1000:.1f} km)."
    )
}

# Ensure results directory exists
out_dir = Path(__file__).resolve().parent.parent.parent / "results"
out_dir.mkdir(exist_ok=True)
out_file = out_dir / "step_045_eccentric_orbit_bound.json"

with open(out_file, "w") as f:
    json.dump(results, f, indent=4)

print(f"Results saved to {out_file}")
print(json.dumps(results, indent=2))
