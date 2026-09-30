#!/usr/bin/env python3
"""
Step 026: Stable Model Comparison with Regularized Likelihood

Four-tier model comparison framework addressing reviewer feedback on
Bayes factor overstatement and parameter-count transparency.

Model Tiers:
------------
1. Null (M_0): 0 parameters, predicts Δv = 0 for all flybys.

2. Anderson Empirical (M_A): 2 parameters. Fitted amplitude on trajectory
   asymmetry (cos δ_in - cos δ_out) plus offset. This is the empirical
   latitude/asymmetry formula from Anderson et al. (2008) simplified to
   asymmetry-only because perigee latitude data are not catalogued.

3. TEP Restricted (M_T^res): 1 parameter. Universal coupling β is the sole
   fitted parameter. All other TEP quantities are pre-specified:
   - λ_TEP ≈ 4200 km (from GNSS atomic clock correlations, Step 016)
   - S_⊕ ≈ 0.35 (from UCD saturation first-principles, Step 010)
   - v_trans ≈ 16.8 km/s (from TEP field equations)
   - Geometry (altitude, velocity, declinations) from JPL Horizons ephemerides

4. TEP Flexible (M_T^flex): 3 parameters. β (conformal coupling),
   b_disf (disformal amplitude), and offset (residual modulation capturing
   plasma screening, OD absorption, and unmodeled systematics).

The restricted model is the scientifically important one because all
parameters except β are pre-specified from independent measurements or
first-principles derivations. The Bayes factor reported for TEP vs Null
refers strictly to the restricted model.

Fixes from prior version:
-------------------------
- Replace the overfitted 4-parameter independent-β_i "Empirical" model with
  the Anderson asymmetry formula.
- Rename "TEP" to "TEP restricted" and add "TEP flexible" as a separate tier.
- Explicitly document pre-specified vs fitted parameters.
- Use log-likelihood and AICc/BIC corrections for stability.

Likelihood / σ_sys (fractional misspecification floor — headline):
------------------------------------------------------------------
**Headline** reported log-likelihoods use a per-flyby fractional
misspecification floor,

    sigma_sys_i = CV_beta * |y_hat_i|,

where y_hat_i is each model's own fitted prediction and
CV_beta = tau / beta_bar is the coefficient of variation of the measured
between-flyby amplitude heterogeneity (DerSimonian–Laird tau on the
per-flyby scale factors, computed in _amplitude_heterogeneity_cv). Each
signal-bearing model is iterated to a fixed point: prediction -> floor ->
refit -> floor, until convergence.

Rationale: the earlier constant geometry-spread floor
(sigma_geom = std of reference-scale predictions at beta_ref = 1e-4) was
anchor-dependent — rescaling the arbitrary normalization beta_ref rescales
the floor and can reorder the model comparison without changing the data.
The fractional floor is dimensionless and invariant under beta_ref; it
encodes that model misspecification scales with the size of the claimed
signal, with the coefficient measured from the catalogue itself.

The Null model predicts zero signal and therefore carries no floor: its
published per-flyby uncertainties are taken at face value, so large
detections penalize it directly.

Sensitivity analyses retained (labelled, non-headline):
- ``likelihood_sensitivity_published_uncertainties_only``: sigma_sys = 0.
  The tiny published uncertainties (~0.01–0.05 mm/s) are inconsistent with
  the ~1–10 mm/s residuals of the single-beta restricted scaling model, so
  this block produces astronomically large BIC values; it is a consistency
  check showing the data taken at face value, not the primary claim.
- ``sensitivity_constant_floor_reference_spread``: the former headline
  constant sigma_geom block, retained to document the anchor-dependent
  convention it superseded.
- ``sensitivity_constant_floor_sweep``: BIC ordering across a grid of
  constant floors, showing where the TEP restricted / Anderson ordering
  crosses zero.
"""

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.flyby_ensemble import (
    ENSEMBLE_GATE_POLICY,
    flyby_ensemble_exclusion_reason,
    flyby_sign_product,
    flyby_snr,
    strict_sign_gate_from_config,
)
from scripts.utils.physics import CHARACTERISTIC_SUPPRESSION
from scripts.utils.step_logger import StepLogger


def convert_to_native_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_to_native_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native_types(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return convert_to_native_types(obj.tolist())
    else:
        return obj


class StableModelComparison:
    """
    Stable model comparison with regularized likelihood.
    """
    
    def __init__(self, logger):
        self.logger = logger

    def _archival_catalog_by_mission(self) -> dict:
        """Map mission_name → flyby dict from Step 003 (provenance for exclusions)."""
        path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
        if not path.is_file():
            return {}
        with open(path, "r", encoding="utf-8") as f:
            catalog = json.load(f)
        out = {}
        for fb in catalog.get("flybys", []):
            name = fb.get("mission_name")
            if name:
                out[str(name)] = fb
        return out

    def load_data(self):
        """
        Load TEP predictions and split catalog entries by Step 008 ensemble gates.

        Returns:
            dict with ``eligible_flybys`` (list for likelihood) and
            ``excluded_from_ensemble`` (audit trail; excluded rows are not dropped
            from the project, only from the gated ensemble).
        """
        pred_file = PROJECT_ROOT / 'results' / 'step007_tep_predictions.json'
        
        if not pred_file.exists():
            self.logger.error(f"TEP predictions not found: {pred_file}")
            return None
        
        with open(pred_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        archival = self._archival_catalog_by_mission()

        flybys = []
        excluded = []
        for name, pred in data['predictions'].items():
            exclusion_reason = flyby_ensemble_exclusion_reason(pred)
            observed = pred.get('observed', {})
            dv_obs = observed.get('dv_obs_mm_s')
            dv_unc = observed.get('sigma_mm_s')
            dv_tep = pred.get('tep_predictions', {}).get('dv_tep_mm_s')
            cat = archival.get(name, {})
            row = {
                'name': name,
                'exclusion_reason': exclusion_reason,
                'dv_obs_mm_s': dv_obs,
                'sigma_mm_s': dv_unc,
                'snr': flyby_snr(dv_obs, dv_unc),
                'dv_tep_mm_s_reference': dv_tep,
                'sign_product': flyby_sign_product(dv_tep, dv_obs),
                'archival_anomaly_reference': cat.get('anomaly_reference'),
                'archival_anomaly_reference_doi': cat.get('anomaly_reference_doi'),
                'archival_usability_notes': cat.get('usability_notes'),
                'archival_detection_significance': cat.get('detection_significance'),
            }
            if exclusion_reason:
                excluded.append(row)
                self.logger.info(
                    f"Skipping {name} for model comparison: {exclusion_reason}"
                )
                continue

            geometry = pred.get('geometry', {})
            tep_pred = pred.get('tep_predictions', {})
            snr = abs(dv_obs) / dv_unc if dv_unc and dv_unc > 0 else 0
            
            flybys.append({
                'name': name,
                'dv_obs': dv_obs,
                'dv_unc': dv_unc,
                'dv_pred_base': tep_pred.get('dv_tep_mm_s', 0.0),
                'dv_grad_base': tep_pred.get('dv_grad_mm_s', 0.0),
                'dv_disf_base': tep_pred.get('dv_disf_mm_s', 0.0),
                'cos_asymmetry': geometry.get('cos_dec_asymmetry', 0.0),
                'snr': snr,
                'archival_anomaly_reference': cat.get('anomaly_reference'),
            })
        
        return {'eligible_flybys': flybys, 'excluded_from_ensemble': excluded}

    def load_full_catalog(self):
        """
        Load ALL flybys with published observations from Step 007, bypassing
        the Step 008 ensemble gates. This expands the model comparison from
        the gated n=6 detections to the full catalog of n=9 flybys with
        explicit TEP predictions and published (or tracking-precision-fallback)
        uncertainties.

        Includes both primary detections and null-result bounds, providing a
        more conservative test because null flybys constrain the models more
        stringently (the Null model predicts zero, which is the observed value).
        """
        pred_file = PROJECT_ROOT / 'results' / 'step007_tep_predictions.json'

        if not pred_file.exists():
            self.logger.error(f"TEP predictions not found: {pred_file}")
            return None

        with open(pred_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        flybys = []
        skipped = []
        for name, pred in data['predictions'].items():
            observed = pred.get('observed', {})
            dv_obs = observed.get('dv_obs_mm_s')
            dv_unc = observed.get('sigma_mm_s')

            # Require both an observed value and an uncertainty
            if dv_obs is None or dv_unc is None or dv_unc <= 0:
                skipped.append({
                    'name': name,
                    'reason': 'missing_observed_or_uncertainty',
                    'dv_obs': dv_obs,
                    'sigma': dv_unc,
                })
                continue

            geometry = pred.get('geometry', {})
            tep_pred = pred.get('tep_predictions', {})
            snr = abs(dv_obs) / dv_unc if dv_unc and dv_unc > 0 else 0

            flybys.append({
                'name': name,
                'dv_obs': dv_obs,
                'dv_unc': dv_unc,
                'dv_pred_base': tep_pred.get('dv_tep_mm_s', 0.0),
                'dv_grad_base': tep_pred.get('dv_grad_mm_s', 0.0),
                'dv_disf_base': tep_pred.get('dv_disf_mm_s', 0.0),
                'cos_asymmetry': geometry.get('cos_dec_asymmetry', 0.0),
                'snr': snr,
            })

        return {'full_catalog_flybys': flybys, 'skipped': skipped}

    def _gauss_loglike(self, residual, sigma_total):
        """Single Gaussian log-likelihood term."""
        return -0.5 * ((residual / sigma_total)**2) - 0.5 * np.log(2 * np.pi * sigma_total**2)

    @staticmethod
    def _sigma_sys_vec(sigma_sys, n):
        """Broadcast a scalar systematic floor or validate a per-flyby vector."""
        arr = np.asarray(sigma_sys, dtype=float)
        if arr.ndim == 0:
            return np.full(n, float(arr))
        if arr.shape != (n,):
            raise ValueError(f"sigma_sys vector must have length {n}, got {arr.shape}")
        return arr

    def log_likelihood_anderson(self, flybys, A, B, sigma_sys=0.0):
        """
        Anderson empirical model: Δv = A * cos_asymmetry + B.
        
        This captures the core empirical finding of Anderson et al. (2008)
        that the anomaly correlates with trajectory asymmetry. Perigee
        latitude is not included because it is not catalogued.
        """
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        log_like = 0.0
        for i, fb in enumerate(flybys):
            dv_pred = A * fb['cos_asymmetry'] + B
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_extra[i]**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_anderson(self, flybys, sigma_sys=0.0):
        """Least-squares fit of Anderson empirical model to data."""
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        X = np.array([[fb['cos_asymmetry'], 1.0] for fb in flybys])
        y = np.array([fb['dv_obs'] for fb in flybys])
        # Weighted least squares with measurement uncertainty
        w = np.array([1.0 / (fb['dv_unc']**2 + sigma_extra[i]**2) for i, fb in enumerate(flybys)])
        W = np.diag(w)
        beta_hat = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
        A_fit, B_fit = beta_hat[0], beta_hat[1]
        log_like = self.log_likelihood_anderson(flybys, A_fit, B_fit, sigma_sys)
        return A_fit, B_fit, log_like

    def log_likelihood_tep_restricted(self, flybys, beta, sigma_sys=0.0):
        """
        TEP Restricted: 1 fitted parameter (β).
        
        All other quantities are pre-specified:
        - λ_TEP ≈ 4000 km (GNSS Step 016)
        - S_⊕ ≈ 0.35 (UCD Step 010)
        - v_trans ≈ 16.8 km/s (TEP field equations)
        - Geometry from JPL Horizons
        
        The TEP velocity shift is linear in β under the canonical amplitude
        sector (Paper 0): the excursion normalization is β-independent, so
            dv_tep ∝ β
        Therefore predictions at arbitrary β are scaled from the reference (β₀=1e-4):
            dv_pred(β) = dv_pred_base * (β / 1e-4)
        """
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        log_like = 0.0
        for i, fb in enumerate(flybys):
            scale = beta / 1e-4
            dv_pred = fb['dv_pred_base'] * scale
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_extra[i]**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_tep_restricted(self, flybys, sigma_sys=0.0):
        """
        Weighted least-squares fit of β for TEP restricted with linear scaling.

        We solve for the optimal scaling factor x = β / 1e-4 directly:
        chi² = Σ (obs - pred_base * x)² / σ²
        d(chi²)/dx = 0  →  x = Σ(pred_base * obs / σ²) / Σ(pred_base² / σ²)
        β_fit = 1e-4 * x

        If x ≤ 0, no positive β reproduces the weighted sign pattern of the
        ensemble; ``beta_fitted`` is None and log-likelihood is evaluated at a
        small positive β for reporting only (near-null TEP), with metadata
        flagging infeasibility of the interior optimum.
        """
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        num = sum(fb['dv_obs'] * fb['dv_pred_base'] / (fb['dv_unc']**2 + sigma_extra[i]**2) for i, fb in enumerate(flybys))
        den = sum(fb['dv_pred_base']**2 / (fb['dv_unc']**2 + sigma_extra[i]**2) for i, fb in enumerate(flybys))
        x = (num / den) if den > 0 else 1.0
        positive_feasible = x > 0
        if positive_feasible:
            beta_fit = 1e-4 * x
            beta_like = beta_fit
        else:
            beta_fit = None
            beta_like = 1e-12
        log_like = self.log_likelihood_tep_restricted(flybys, beta_like, sigma_sys)
        meta = {
            'optimal_scale_factor_x': float(x),
            'positive_beta_scale_feasible': bool(positive_feasible),
            'log_likelihood_evaluated_at_beta': float(beta_like),
        }
        return beta_fit, log_like, meta

    def log_likelihood_tep_flexible(self, flybys, beta, b_disf, offset, sigma_sys=0.0):
        """
        TEP Flexible: 3 fitted parameters (β, b_disf, offset).
        
        Both gradient and disformal components scale linearly in β:
            dv_grad(β) = dv_grad_base * (β / 1e-4)
            dv_disf(β) = dv_disf_base * (β / 1e-4)

        The prediction is:
            Δv = (β / 1e-4) * (dv_grad_base + b_disf * dv_disf_base) + offset
        
        b_disf allows the disformal amplitude to vary independently (as a ratio),
        and offset captures residual modulation (plasma, OD, etc.).
        
        Note: This is a non-linear model in β. For fitting, we use an iterative
        approach: first fit the linear combination at a reference β, then optimize β.
        """
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        log_like = 0.0
        scale = beta / 1e-4
        for i, fb in enumerate(flybys):
            dv_pred = scale * (fb['dv_grad_base'] + b_disf * fb['dv_disf_base']) + offset
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_extra[i]**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_tep_flexible(self, flybys, sigma_sys=0.0):
        """
        Fit TEP flexible with linear (canonical) beta scaling.
        
        Strategy: For a given β, the model is linear in (b_disf, offset) with
        design matrix columns (scale*dv_disf_base, 1.0). We iterate over β to
        find the global optimum.
        """
        from scipy.optimize import minimize_scalar
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))

        def _neg_loglike_at_beta(log_beta):
            beta = np.exp(log_beta)
            scale = beta / 1e-4
            # Linear least squares for b_disf and offset at fixed beta
            X = np.array([[scale * fb['dv_disf_base'], 1.0] for fb in flybys])
            y = np.array([fb['dv_obs'] - scale * fb['dv_grad_base'] for fb in flybys])
            w = np.array([1.0 / (fb['dv_unc']**2 + sigma_extra[i]**2) for i, fb in enumerate(flybys)])
            W = np.diag(w)
            
            try:
                beta_hat = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
                b_disf_fit, offset_fit = beta_hat[0], beta_hat[1]
            except (np.linalg.LinAlgError, ValueError):
                return 1e10
            
            log_like = self.log_likelihood_tep_flexible(flybys, beta, b_disf_fit, offset_fit, sigma_sys)
            return -log_like
        
        # Optimize log_beta in a broad range around the step008 value
        result = minimize_scalar(_neg_loglike_at_beta, bounds=(-12, -6), method='bounded')
        beta_fit = np.exp(result.x)
        
        # Refit b_disf and offset at optimal beta
        scale = beta_fit / 1e-4
        X = np.array([[scale * fb['dv_disf_base'], 1.0] for fb in flybys])
        y = np.array([fb['dv_obs'] - scale * fb['dv_grad_base'] for fb in flybys])
        w = np.array([1.0 / (fb['dv_unc']**2 + sigma_extra[i]**2) for i, fb in enumerate(flybys)])
        W = np.diag(w)
        beta_hat = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
        b_disf_fit, offset_fit = beta_hat[0], beta_hat[1]
        
        log_like = self.log_likelihood_tep_flexible(flybys, beta_fit, b_disf_fit, offset_fit, sigma_sys)
        return beta_fit, b_disf_fit, offset_fit, log_like

    def log_likelihood_null(self, flybys, sigma_sys=0.0):
        """
        Null model: 0 parameters, predicts Δv = 0.
        """
        sigma_extra = self._sigma_sys_vec(sigma_sys, len(flybys))
        log_like = 0.0
        for i, fb in enumerate(flybys):
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_extra[i]**2)
            residual = fb['dv_obs']
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like
    
    def compute_aic_bic(self, log_like, n_params, n_data):
        """
        Compute AIC, AICc, and BIC with proper corrections.
        
        AIC = -2*log(L) + 2*k
        AICc = AIC + 2*k*(k+1)/(n-k-1)  (small sample correction)
        BIC = -2*log(L) + k*ln(n)
        """
        aic = -2 * log_like + 2 * n_params
        
        # AICc correction for small samples
        if n_data - n_params - 1 > 0:
            aic_correction = 2 * n_params * (n_params + 1) / (n_data - n_params - 1)
            aicc = aic + aic_correction
        else:
            aicc = float('inf')  # Not enough data
        
        bic = -2 * log_like + n_params * np.log(n_data)
        
        return aic, aicc, bic
    
    # Threshold above which the BIC approximation is not considered valid for
    # small samples (n < 10).  log10(BF) > 100 is already astronomically large
    # and signals the approximation is being driven by formal uncertainties
    # rather than genuine information-theoretic compression.
    LOG10_BF_VALIDITY_THRESHOLD = 100.0

    def bayes_factor_approx(self, log_like_1, log_like_2, k1, k2, n):
        """
        Approximate Bayes factor using BIC approximation.
        
        BF ≈ exp((BIC2 - BIC1)/2)
        
        This is a large-sample approximation and breaks down for very small
        samples (n < 10) or extreme signal-to-noise ratios.  When log10(BF)
        exceeds ~100 the number is driven by formal uncertainties and should
        not be reported as a literal probability ratio.
        """
        bic1 = -2 * log_like_1 + k1 * np.log(n)
        bic2 = -2 * log_like_2 + k2 * np.log(n)
        
        delta_bic = bic2 - bic1
        log_bf = delta_bic / 2
        log10_bf = float(log_bf / math.log(10.0))

        # float64 overflows at exp(709.78...).
        if log_bf > 709.0:
            bf = float("inf")
        elif log_bf < -745.0:
            bf = 0.0
        else:
            bf = float(np.exp(log_bf))

        # Flag when the BIC approximation is outside its domain of validity.
        # For n < 10 the large-sample Laplace approximation is unreliable,
        # and log10(BF) > 100 is a clear signal of formal-uncertainty dominance.
        approximation_valid = (
            n >= 10 and abs(log10_bf) <= self.LOG10_BF_VALIDITY_THRESHOLD
        )

        return bf, float(delta_bic), log10_bf, approximation_valid

    def _evaluate_four_tiers(self, flybys, sigma_sys: float, *, log_details: bool = True):
        """
        Run all four likelihood tiers at a fixed extra Gaussian σ_sys added in
        quadrature to each flyby's published uncertainty.
        """
        n_data = len(flybys)
        log_like_null = self.log_likelihood_null(flybys, sigma_sys)
        aic_null, aicc_null, bic_null = self.compute_aic_bic(log_like_null, 0, n_data)

        A_fit, B_fit, log_like_anderson = self.fit_anderson(flybys, sigma_sys)
        aic_anderson, aicc_anderson, bic_anderson = self.compute_aic_bic(log_like_anderson, 2, n_data)

        beta_fit, log_like_tep_res, tep_res_meta = self.fit_tep_restricted(flybys, sigma_sys)
        aic_tep_res, aicc_tep_res, bic_tep_res = self.compute_aic_bic(log_like_tep_res, 1, n_data)

        beta_flex, b_disf_flex, offset_flex, log_like_tep_flex = self.fit_tep_flexible(flybys, sigma_sys)
        aic_tep_flex, aicc_tep_flex, bic_tep_flex = self.compute_aic_bic(log_like_tep_flex, 3, n_data)

        if log_details:
            self.logger.info("")
            self.logger.info("Fitted parameters:")
            self.logger.info(f"  Anderson: A={A_fit:.3f}, B={B_fit:.3f}")
            br = f"{beta_fit:.2e}" if beta_fit is not None else "None (no positive-β interior optimum)"
            self.logger.info(f"  TEP restricted: beta={br}")
            self.logger.info(
                f"  TEP flexible: beta={beta_flex:.2e}, b_disf={b_disf_flex:.2e}, offset={offset_flex:.3f}"
            )

            self.logger.info("")
            self.logger.info("Log-likelihoods:")
            self.logger.info(f"  Null:        {log_like_null:.2f}")
            self.logger.info(f"  Anderson:    {log_like_anderson:.2f}")
            self.logger.info(f"  TEP res:     {log_like_tep_res:.2f}")
            self.logger.info(f"  TEP flex:    {log_like_tep_flex:.2f}")

            self.logger.info("")
            self.logger.info("Information Criteria:")
            self.logger.info(f"  Null:        AIC={aic_null:.1f}, BIC={bic_null:.1f}")
            self.logger.info(f"  Anderson:    AIC={aic_anderson:.1f}, BIC={bic_anderson:.1f}")
            self.logger.info(f"  TEP res:     AIC={aic_tep_res:.1f}, BIC={bic_tep_res:.1f}")
            self.logger.info(f"  TEP flex:    AIC={aic_tep_flex:.1f}, BIC={bic_tep_flex:.1f}")

        # Bayes factors (all relative to Null)
        bf_anderson_null, dbic_anderson_null, log10_bf_anderson_null, valid_anderson_null = self.bayes_factor_approx(
            log_like_anderson, log_like_null, 2, 0, n_data)
        bf_tep_res_null, dbic_tep_res_null, log10_bf_tep_res_null, valid_tep_res_null = self.bayes_factor_approx(
            log_like_tep_res, log_like_null, 1, 0, n_data)
        bf_tep_flex_null, dbic_tep_flex_null, log10_bf_tep_flex_null, valid_tep_flex_null = self.bayes_factor_approx(
            log_like_tep_flex, log_like_null, 3, 0, n_data)

        # Anderson vs TEP restricted
        bf_tep_res_anderson, dbic_tep_res_anderson, log10_bf_tep_res_anderson, valid_tep_res_anderson = self.bayes_factor_approx(
            log_like_tep_res, log_like_anderson, 1, 2, n_data)

        if log_details:
            self.logger.info("")
            self.logger.info("Bayes Factors (via BIC; log10 shown):")
            # Warn when the BIC approximation is outside its validity regime.
            any_invalid = not (
                valid_anderson_null and valid_tep_res_null
                and valid_tep_flex_null and valid_tep_res_anderson
            )
            if any_invalid:
                self.logger.info(
                    f"  WARNING: n={n_data} is too small for the BIC large-sample approximation."
                )
                self.logger.info(
                    f"  log10(BF) values exceeding {self.LOG10_BF_VALIDITY_THRESHOLD} are driven by formal"
                )
                self.logger.info(
                    "  uncertainties and should not be reported as literal probability ratios."
                )
            self.logger.info(
                f"  Anderson vs Null:       log10(B)={log10_bf_anderson_null:.2f} (ΔBIC={dbic_anderson_null:.1f})"
            )
            self.logger.info(
                f"  TEP restricted vs Null: log10(B)={log10_bf_tep_res_null:.2f} (ΔBIC={dbic_tep_res_null:.1f})"
            )
            self.logger.info(
                f"  TEP flexible vs Null:   log10(B)={log10_bf_tep_flex_null:.2f} (ΔBIC={dbic_tep_flex_null:.1f})"
            )
            self.logger.info(
                f"  TEP res vs Anderson:    log10(B)={log10_bf_tep_res_anderson:.2f} (ΔBIC={dbic_tep_res_anderson:.1f})"
            )
        
        # Model selection based on BIC
        bic_values = {
            'Null': bic_null,
            'Anderson': bic_anderson,
            'TEP_restricted': bic_tep_res,
            'TEP_flexible': bic_tep_flex
        }
        best_model_bic = min(bic_values, key=bic_values.get)
        
        # Akaike weights (TEP restricted vs Null vs Anderson)
        aicc_comp = {'TEP_restricted': aicc_tep_res, 'Null': aicc_null, 'Anderson': aicc_anderson}
        finite_aicc = {k: v for k, v in aicc_comp.items() if np.isfinite(v)}
        if finite_aicc:
            min_aicc = min(finite_aicc.values())
            delta_aicc = {k: v - min_aicc for k, v in finite_aicc.items()}
            sum_exp = sum(np.exp(-0.5 * da) for da in delta_aicc.values())
            akaike_weights = {k: np.exp(-0.5 * da) / sum_exp for k, da in delta_aicc.items()}
        else:
            akaike_weights = {'TEP_restricted': 1.0/3, 'Null': 1.0/3, 'Anderson': 1.0/3}
        
        results = {
            'n_data': n_data,
            'extra_sigma_sys_quadrature_mm_s': float(sigma_sys),
            'log_likelihoods': {
                'Null': float(log_like_null),
                'Anderson': float(log_like_anderson),
                'TEP_restricted': float(log_like_tep_res),
                'TEP_flexible': float(log_like_tep_flex)
            },
            'fitted_parameters': {
                'Anderson': {'A': float(A_fit), 'B': float(B_fit)},
                'TEP_restricted': {
                    'beta': float(beta_fit) if beta_fit is not None else None,
                    'fit_diagnostics': tep_res_meta,
                },
                'TEP_flexible': {
                    'beta': float(beta_flex),
                    'b_disf': float(b_disf_flex),
                    'offset': float(offset_flex)
                }
            },
            'information_criteria': {
                'Null': {'AIC': float(aic_null), 'AICc': float(aicc_null), 'BIC': float(bic_null)},
                'Anderson': {'AIC': float(aic_anderson), 'AICc': float(aicc_anderson), 'BIC': float(bic_anderson)},
                'TEP_restricted': {'AIC': float(aic_tep_res), 'AICc': float(aicc_tep_res), 'BIC': float(bic_tep_res)},
                'TEP_flexible': {'AIC': float(aic_tep_flex), 'AICc': float(aicc_tep_flex), 'BIC': float(bic_tep_flex)}
            },
            'bayes_factors': {
                'Anderson_vs_Null': bf_anderson_null,
                'TEP_restricted_vs_Null': bf_tep_res_null,
                'TEP_flexible_vs_Null': bf_tep_flex_null,
                'TEP_restricted_vs_Anderson': bf_tep_res_anderson,
                'delta_BIC_Anderson_vs_Null': dbic_anderson_null,
                'delta_BIC_TEP_restricted_vs_Null': dbic_tep_res_null,
                'delta_BIC_TEP_flexible_vs_Null': dbic_tep_flex_null,
                'delta_BIC_TEP_restricted_vs_Anderson': dbic_tep_res_anderson,
                'log10_BF_Anderson_vs_Null': log10_bf_anderson_null,
                'log10_BF_TEP_restricted_vs_Null': log10_bf_tep_res_null,
                'log10_BF_TEP_flexible_vs_Null': log10_bf_tep_flex_null,
                'log10_BF_TEP_restricted_vs_Anderson': log10_bf_tep_res_anderson,
                'bic_approximation_valid': {
                    'Anderson_vs_Null': valid_anderson_null,
                    'TEP_restricted_vs_Null': valid_tep_res_null,
                    'TEP_flexible_vs_Null': valid_tep_flex_null,
                    'TEP_restricted_vs_Anderson': valid_tep_res_anderson,
                },
                'bic_approximation_note': (
                    "The BIC approximation BF ≈ exp(ΔBIC/2) is a large-sample result. "
                    "For n < 10 and extreme S/N it is unreliable.  When log10(BF) > 100 "
                    "the value is driven by formal uncertainties and should not be "
                    "reported as a literal probability ratio."
                ),
            },
            'model_selection': {
                'best_model_BIC': best_model_bic,
                'akaike_weights': akaike_weights,
                'interpretation': f"BIC selects {best_model_bic} model"
            },
            'pre_specified_parameters': {
                'lambda_TEP_km': 4200,
                'lambda_TEP_source': 'GNSS atomic clock correlations (Step 016)',
                'S_earth': round(CHARACTERISTIC_SUPPRESSION, 3),
                'S_earth_source': 'UCD saturation first-principles (Step 010)',
                'v_trans_km_s': 16.8,
                'v_trans_source': 'TEP field equations',
                'geometry_source': 'JPL Horizons ephemerides'
            }
        }

        return results

    def stable_model_comparison(self, flybys):
        """
        Headline comparison uses the fractional misspecification floor
        (sigma_sys_i = CV_beta * |y_hat_i| per flyby, iterated to a fixed
        point), which is invariant under the arbitrary beta_ref anchor.
        The former constant geometry-spread floor and a
        published-uncertainties-only block (sigma_sys = 0) are retained as
        labelled sensitivity analyses, not the primary evidence claim.
        """
        self.logger.section("FOUR-TIER MODEL COMPARISON")

        n_data = len(flybys)
        self.logger.info(f"Using {n_data} flybys passing Step 008 ensemble gates")

        pred_vec = np.array([fb["dv_pred_base"] for fb in flybys], dtype=float)
        sigma_geom = float(np.std(pred_vec, ddof=1)) if len(pred_vec) > 1 else 0.0

        self.logger.info(
            f"Anchor-dependent reference spread σ_geom = {sigma_geom:.3f} mm/s "
            "(constant floor, sensitivity only)."
        )
        self.logger.info(
            "The constant floor is retained only as a labelled sensitivity "
            "because it rescales with the arbitrary beta_ref anchor."
        )

        frac = self.fractional_floor_model_comparison(flybys)
        if frac is not None:
            primary = frac
            self.logger.info(
                f"Headline systematic model: fractional misspecification floor "
                f"(cv_amp = {frac['amplitude_heterogeneity']['cv_amp']:.3f})."
            )
        else:
            primary = self._evaluate_four_tiers(flybys, sigma_geom, log_details=True)
            primary["systematic_uncertainty_mm_s"] = sigma_geom
            primary["systematic_uncertainty_model"] = "geometry_prediction_spread_std_ddof1"
        primary["geometry_prediction_spread_mm_s"] = sigma_geom
        primary["reference_scale_spread_mm_s_anchor_dependent"] = sigma_geom

        # Sensitivity: published uncertainties only (sigma_sys = 0)
        # These produce extremely large BIC values because tiny published
        # uncertainties are inconsistent with single-β restricted scaling residuals.
        sens = self._evaluate_four_tiers(flybys, 0.0, log_details=False)
        primary["likelihood_sensitivity_published_uncertainties_only"] = {
            "extra_sigma_sys_quadrature_mm_s": 0.0,
            "note": (
                "Published per-flyby uncertainties without any additional systematic term. "
                "These produce extremely large BIC values because the tiny published "
                "uncertainties (~0.01–0.05 mm/s) are inconsistent with the ~1–10 mm/s "
                "residuals of the single-β restricted scaling model. Reported for transparency and "
                "consistency checking only; the fractional-floor comparison is the "
                "scientifically meaningful one."
            ),
            "log_likelihoods": sens["log_likelihoods"],
            "information_criteria": sens["information_criteria"],
            "bayes_factors": sens["bayes_factors"],
        }

        return primary

    def ensemble_composition_robustness(self, eligible_flybys, excluded_catalog, primary_results):
        """
        Robustness test: evaluate model comparison under alternative ensemble
        compositions to address reviewer concerns about sample-size sensitivity.

        Specifically tests the historical n=4 composition (Cassini re-added
        despite sign mismatch under the pre-retranscription catalogue) vs the
        headline gated ensemble when ``strict_sign_gate`` is true in
        ``config/pipeline_config.json``. When Cassini is already in ``eligible_flybys``,
        this block is skipped as degenerate.
        """
        self.logger.section("ENSEMBLE COMPOSITION ROBUSTNESS")
        self.logger.info(
            "Testing whether model-comparison conclusions are sensitive to "
            "Cassini exclusion vs inclusion (historical n vs n+1 composition)."
        )

        names = {fb["name"] for fb in eligible_flybys}
        if "Cassini" in names:
            self.logger.info(
                "Cassini already in eligible ensemble (strict_sign_gate=false); "
                "skipping separate n=4 re-merge."
            )
            return {
                "skipped": True,
                "reason": "cassini_already_in_primary_ensemble",
                "note": (
                    "With strict_sign_gate disabled, Step 026's primary gated row set "
                    "already includes Cassini; the historical n=4 robustness merge is redundant."
                ),
            }

        # Build n=4 list by adding Cassini back from excluded catalog
        cassini_row = None
        for row in excluded_catalog:
            if row.get("name") == "Cassini":
                cassini_row = row
                break

        if cassini_row is None:
            self.logger.info("Cassini not found in excluded catalog; skipping n=4 test.")
            return None

        pred_file = PROJECT_ROOT / "results" / "step007_tep_predictions.json"
        cassini_pred = {}
        if pred_file.is_file():
            with open(pred_file, "r", encoding="utf-8") as f:
                cassini_pred = json.load(f).get("predictions", {}).get("Cassini", {})

        tep = cassini_pred.get("tep_predictions", {}) if isinstance(cassini_pred, dict) else {}
        dv_grad = float(tep.get("dv_grad_mm_s", 0.0) or 0.0)
        dv_disf = float(tep.get("dv_disf_mm_s", 0.0) or 0.0)
        cos_asym = float(cassini_pred.get("geometry", {}).get("cos_dec_asymmetry", 0.0) or 0.0)

        # Reconstruct Cassini as a flyby dict matching the eligible format
        cassini_flyby = {
            "name": cassini_row["name"],
            "dv_obs": cassini_row["dv_obs_mm_s"],
            "dv_unc": cassini_row["sigma_mm_s"],
            "dv_pred_base": cassini_row.get("dv_tep_mm_s_reference", 0.0),
            "dv_grad_base": dv_grad,
            "dv_disf_base": dv_disf,
            "cos_asymmetry": cos_asym,
            "snr": cassini_row.get("snr", 0.0),
        }

        n4_flybys = eligible_flybys + [cassini_flyby]

        # Use the same fractional misspecification floor as the headline comparison
        n4_result = self.fractional_floor_model_comparison(n4_flybys)
        if n4_result is None:
            self.logger.info("Fractional-floor comparison unavailable for n=4; skipping.")
            return None

        self.logger.info(f"Headline n={len(eligible_flybys)} vs robustness n={len(n4_flybys)}")
        self.logger.info(
            f"  TEP restricted vs Null log10(BF): n={len(eligible_flybys)} = "
            f"{primary_results['bayes_factors']['log10_BF_TEP_restricted_vs_Null']:.2f}; "
            f"n={len(n4_flybys)} = {n4_result['bayes_factors']['log10_BF_TEP_restricted_vs_Null']:.2f}"
        )
        self.logger.info(
            f"  Best model BIC: n={len(eligible_flybys)} = {primary_results['model_selection']['best_model_BIC']}; "
            f"n={len(n4_flybys)} = {n4_result['model_selection']['best_model_BIC']}"
        )

        return {
            "note": (
                "Model comparison evaluated on n=4 (Cassini included) to test "
                "sensitivity of conclusions to the sign-mismatch exclusion. "
                "Uses the same fractional misspecification floor as the headline comparison."
            ),
            "n_data": len(n4_flybys),
            "cassini_included": True,
            "systematic_uncertainty_model": n4_result.get("systematic_uncertainty_model"),
            "log_likelihoods": n4_result["log_likelihoods"],
            "information_criteria": n4_result["information_criteria"],
            "bayes_factors": n4_result["bayes_factors"],
            "model_selection": n4_result["model_selection"],
            "fitted_parameters": n4_result["fitted_parameters"],
        }

    def full_catalog_model_comparison(self, full_catalog_flybys):
        """
        Breakthrough test: four-tier model comparison on the FULL catalog of
        all flybys with published observations and explicit TEP predictions.

        Expands the sample from the gated n=6 to n=9, including:
        - 6 published anomalies (NEAR, Galileo 1990, Rosetta 2005, Cassini, Galileo 1992, MESSENGER)
        - 3 null-result bounds (Juno, Rosetta 2007, Rosetta 2009)

        This is a more conservative test than the gated ensemble because:
        1. The BIC large-sample approximation is more reliable at n=9.
        2. Null flybys provide constraints that penalize overfitting.
        3. Cassini's sign mismatch must be explained by the TEP flexible model.

        The headline systematic floor is the fractional-misspecification
        model (``fractional_floor_model_comparison``): each signal-bearing
        model carries a per-flyby systematic proportional to its own fitted
        amplitude, with the proportionality measured from the catalogue's
        own between-flyby amplitude heterogeneity.  The legacy constant
        floor sigma_geom = std(reference-scale predictions) is retained as a
        labelled sensitivity, alongside a constant-floor sweep, because it
        is anchor-dependent (it rescales with beta_ref) and must not be the
        headline.
        """
        self.logger.section("FULL-CATALOG MODEL COMPARISON (n=9)")
        self.logger.info(
            "Testing whether model-preference conclusions hold when the "
            "comparison is forced to explain ALL published flybys, not only "
            "the three sign-gated primary detections."
        )

        n_data = len(full_catalog_flybys)
        self.logger.info(f"Using {n_data} flybys with published observations from Step 007")

        pred_vec = np.array([fb["dv_pred_base"] for fb in full_catalog_flybys], dtype=float)
        sigma_geom = float(np.std(pred_vec, ddof=1)) if len(pred_vec) > 1 else 0.0

        # Primary: anchor-invariant fractional misspecification floor
        frac = self.fractional_floor_model_comparison(full_catalog_flybys)
        if frac is not None:
            result = frac
            self.logger.info(
                f"Headline systematic model: fractional misspecification floor "
                f"(cv_amp = {frac['amplitude_heterogeneity']['cv_amp']:.3f}); "
                "per-flyby floor scales with each model's fitted signal, "
                "invariant under the beta_ref anchor."
            )
        else:
            self.logger.warning(
                "Fractional-floor comparison unavailable; falling back to "
                "constant sigma_geom floor (anchor-dependent)."
            )
            result = self._evaluate_four_tiers(
                full_catalog_flybys, sigma_geom, log_details=True)
            result["systematic_uncertainty_mm_s"] = sigma_geom
            result["systematic_uncertainty_model"] = \
                "geometry_prediction_spread_full_catalog_std_ddof1"

        result["n_data"] = n_data
        result["included_flybys"] = [fb["name"] for fb in full_catalog_flybys]
        result["included_flyby_summary"] = {
            "n_detections": sum(1 for fb in full_catalog_flybys if fb["snr"] >= 2),
            "n_nulls": sum(1 for fb in full_catalog_flybys if fb["snr"] < 2),
        }

        # Sensitivity: legacy anchor-dependent constant floor
        legacy = self._evaluate_four_tiers(
            full_catalog_flybys, sigma_geom, log_details=False)
        result["sensitivity_constant_floor_reference_spread"] = {
            "systematic_uncertainty_mm_s": sigma_geom,
            "note": (
                "Constant floor sigma_sys = std(reference-scale predictions). "
                "Anchor-dependent: it rescales with the arbitrary beta_ref "
                "anchor, so it is retained as a sensitivity diagnostic only."
            ),
            "information_criteria": legacy["information_criteria"],
            "bayes_factors": legacy["bayes_factors"],
            "log_likelihoods": legacy["log_likelihoods"],
        }

        # Sensitivity: constant-floor ordering sweep
        sweep_floors = sorted({
            0.0, sigma_geom, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0
        })
        result["sensitivity_constant_floor_sweep"] = \
            self.constant_floor_sensitivity_sweep(full_catalog_flybys, sweep_floors)

        self.logger.info("")
        self.logger.info("=" * 70)
        self.logger.info("FULL-CATALOG RESULTS SUMMARY")
        self.logger.info(f"  n_data = {n_data}")
        self.logger.info(
            f"  TEP restricted vs Null log10(BF) = "
            f"{result['bayes_factors']['log10_BF_TEP_restricted_vs_Null']:.2f}"
        )
        self.logger.info(
            f"  TEP restricted vs Anderson log10(BF) = "
            f"{result['bayes_factors']['log10_BF_TEP_restricted_vs_Anderson']:.2f}"
        )
        self.logger.info(
            f"  Best model (BIC) = {result['model_selection']['best_model_BIC']}"
        )
        self.logger.info("=" * 70)

        return result

    def _amplitude_heterogeneity_cv(self, flybys):
        """
        Fractional between-flyby scatter of the fitted response amplitude,
        estimated by DerSimonian–Laird on the implied scale factors
        x_i = dv_obs_i / dv_pred_base_i over the amplitude-informative rows
        (|dv_pred_base| > 0).  Returns cv_amp = tau_x / |x_bar|.

        This is a data property — it measures how much the single-amplitude
        restricted model is misspecified per flyby — and is invariant under
        the arbitrary reference coupling beta_ref because it is a ratio of
        fitted-amplitude quantities.
        """
        base = np.array([fb["dv_pred_base"] for fb in flybys], dtype=float)
        obs = np.array([fb["dv_obs"] for fb in flybys], dtype=float)
        unc = np.array([fb["dv_unc"] for fb in flybys], dtype=float)
        inform = np.abs(base) > 1e-8
        if inform.sum() < 2:
            return None
        x_i = obs[inform] / base[inform]
        s_i = unc[inform] / np.abs(base[inform])
        w = 1.0 / s_i**2
        x_bar = float(np.sum(w * x_i) / np.sum(w))
        Q = float(np.sum(w * (x_i - x_bar) ** 2))
        df = int(inform.sum()) - 1
        C = float(np.sum(w) - np.sum(w**2) / np.sum(w))
        tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
        return {
            "cv_amp": float(np.sqrt(tau2) / abs(x_bar)) if x_bar != 0 else None,
            "tau_scale": float(np.sqrt(tau2)),
            "x_bar_scale": x_bar,
            "n_informative": int(inform.sum()),
            "cochran_Q": Q,
        }

    def fractional_floor_model_comparison(self, flybys):
        """
        Anchor-invariant four-tier comparison using a fractional
        misspecification floor: each signal-bearing model carries a per-flyby
        systematic proportional to its own fitted signal amplitude,

            sigma_sys_i = cv_amp * |pred_i(model)|

        so the floor scales with the signal the model asserts rather than
        with the arbitrary reference coupling beta_ref.  cv_amp is measured
        from the catalogue itself via DerSimonian–Laird on the implied
        per-flyby amplitude scale factors.  The Null model predicts zero
        signal, hence zero floor: it is evaluated against the published
        uncertainties at face value (the anomalies are real detections, so
        this is the correct honest penalty).  Each model is iterated to the
        fixed point pred -> floor -> refit.
        """
        het = self._amplitude_heterogeneity_cv(flybys)
        if het is None or het["cv_amp"] is None:
            return None
        cv = het["cv_amp"]
        n_data = len(flybys)
        base = np.array([fb["dv_pred_base"] for fb in flybys], dtype=float)
        asym = np.array([fb["cos_asymmetry"] for fb in flybys], dtype=float)
        grad = np.array([fb["dv_grad_base"] for fb in flybys], dtype=float)
        disf = np.array([fb["dv_disf_base"] for fb in flybys], dtype=float)
        obs = np.array([fb["dv_obs"] for fb in flybys], dtype=float)

        # Null: zero predicted signal -> zero misspecification floor
        log_like_null = self.log_likelihood_null(flybys, 0.0)

        # TEP restricted fixed point
        pred_tr = np.zeros(n_data)
        for _ in range(80):
            floor = cv * np.abs(pred_tr)
            beta_tr, _, _ = self.fit_tep_restricted(flybys, floor)
            scale = (beta_tr / 1e-4) if beta_tr else 1e-8
            pred_new = base * scale
            if np.allclose(pred_new, pred_tr, atol=1e-12):
                pred_tr = pred_new
                break
            pred_tr = pred_new
        floor_tr = cv * np.abs(pred_tr)
        beta_tr, log_like_tr, meta_tr = self.fit_tep_restricted(flybys, floor_tr)

        # Anderson fixed point
        pred_and = np.zeros(n_data)
        for _ in range(80):
            floor = cv * np.abs(pred_and)
            A_fit, B_fit, _ = self.fit_anderson(flybys, floor)
            pred_new = A_fit * asym + B_fit
            if np.allclose(pred_new, pred_and, atol=1e-12):
                pred_and = pred_new
                break
            pred_and = pred_new
        floor_and = cv * np.abs(pred_and)
        A_fit, B_fit, log_like_and = self.fit_anderson(flybys, floor_and)

        # TEP flexible fixed point
        pred_tf = np.zeros(n_data)
        beta_flex, b_disf_flex, offset_flex = 1e-4, 0.0, 0.0
        for _ in range(80):
            floor = cv * np.abs(pred_tf)
            beta_flex, b_disf_flex, offset_flex, _ = self.fit_tep_flexible(
                flybys, floor)
            pred_new = (beta_flex / 1e-4) * (grad + b_disf_flex * disf) + offset_flex
            if np.allclose(pred_new, pred_tf, atol=1e-12):
                pred_tf = pred_new
                break
            pred_tf = pred_new
        floor_tf = cv * np.abs(pred_tf)
        beta_flex, b_disf_flex, offset_flex, log_like_tf = \
            self.fit_tep_flexible(flybys, floor_tf)

        aic_null, aicc_null, bic_null = self.compute_aic_bic(log_like_null, 0, n_data)
        aic_and, aicc_and, bic_and = self.compute_aic_bic(log_like_and, 2, n_data)
        aic_tr, aicc_tr, bic_tr = self.compute_aic_bic(log_like_tr, 1, n_data)
        aic_tf, aicc_tf, bic_tf = self.compute_aic_bic(log_like_tf, 3, n_data)

        bf_tr_null, dbic_tr_null, l10_tr_null, valid_tr_null = \
            self.bayes_factor_approx(log_like_tr, log_like_null, 1, 0, n_data)
        bf_and_null, dbic_and_null, l10_and_null, valid_and_null = \
            self.bayes_factor_approx(log_like_and, log_like_null, 2, 0, n_data)
        bf_tf_null, dbic_tf_null, l10_tf_null, valid_tf_null = \
            self.bayes_factor_approx(log_like_tf, log_like_null, 3, 0, n_data)
        bf_tr_and, dbic_tr_and, l10_tr_and, valid_tr_and = \
            self.bayes_factor_approx(log_like_tr, log_like_and, 1, 2, n_data)
        bf_tf_and, dbic_tf_and, l10_tf_and, valid_tf_and = \
            self.bayes_factor_approx(log_like_tf, log_like_and, 3, 2, n_data)

        bic_values = {
            'Null': bic_null, 'Anderson': bic_and,
            'TEP_restricted': bic_tr, 'TEP_flexible': bic_tf,
        }
        best_model_bic = min(bic_values, key=bic_values.get)

        aicc_comp = {
            'TEP_restricted': aicc_tr, 'Null': aicc_null,
            'Anderson': aicc_and, 'TEP_flexible': aicc_tf,
        }
        finite_aicc = {k: v for k, v in aicc_comp.items() if np.isfinite(v)}
        if finite_aicc:
            min_aicc = min(finite_aicc.values())
            delta_aicc = {k: v - min_aicc for k, v in finite_aicc.items()}
            sum_exp = sum(np.exp(-0.5 * da) for da in delta_aicc.values())
            akaike_weights = {k: np.exp(-0.5 * da) / sum_exp for k, da in delta_aicc.items()}
        else:
            akaike_weights = {k: 0.25 for k in aicc_comp}

        return {
            "systematic_uncertainty_model": (
                "fractional_misspecification_floor_cv_times_abs_prediction"
            ),
            "amplitude_heterogeneity": het,
            "n_data": n_data,
            "per_flyby_floor_mm_s": {
                "TEP_restricted": {
                    fb["name"]: float(v) for fb, v in zip(flybys, floor_tr)
                },
                "Anderson": {
                    fb["name"]: float(v) for fb, v in zip(flybys, floor_and)
                },
                "TEP_flexible": {
                    fb["name"]: float(v) for fb, v in zip(flybys, floor_tf)
                },
            },
            "log_likelihoods": {
                'Null': float(log_like_null),
                'Anderson': float(log_like_and),
                'TEP_restricted': float(log_like_tr),
                'TEP_flexible': float(log_like_tf),
            },
            "fitted_parameters": {
                'Anderson': {'A': float(A_fit), 'B': float(B_fit)},
                'TEP_restricted': {
                    'beta': float(beta_tr) if beta_tr is not None else None,
                    'fit_diagnostics': meta_tr,
                },
                'TEP_flexible': {
                    'beta': float(beta_flex),
                    'b_disf': float(b_disf_flex),
                    'offset': float(offset_flex),
                },
            },
            "information_criteria": {
                'Null': {'AIC': aic_null, 'AICc': aicc_null, 'BIC': bic_null},
                'Anderson': {'AIC': aic_and, 'AICc': aicc_and, 'BIC': bic_and},
                'TEP_restricted': {'AIC': aic_tr, 'AICc': aicc_tr, 'BIC': bic_tr},
                'TEP_flexible': {'AIC': aic_tf, 'AICc': aicc_tf, 'BIC': bic_tf},
            },
            "bayes_factors": {
                'Anderson_vs_Null': bf_and_null,
                'TEP_restricted_vs_Null': bf_tr_null,
                'TEP_flexible_vs_Null': bf_tf_null,
                'TEP_restricted_vs_Anderson': bf_tr_and,
                'TEP_flexible_vs_Anderson': bf_tf_and,
                'delta_BIC_Anderson_vs_Null': dbic_and_null,
                'delta_BIC_TEP_restricted_vs_Null': dbic_tr_null,
                'delta_BIC_TEP_flexible_vs_Null': dbic_tf_null,
                'delta_BIC_TEP_restricted_vs_Anderson': dbic_tr_and,
                'delta_BIC_TEP_flexible_vs_Anderson': dbic_tf_and,
                'log10_BF_Anderson_vs_Null': l10_and_null,
                'log10_BF_TEP_restricted_vs_Null': l10_tr_null,
                'log10_BF_TEP_restricted_vs_Anderson': l10_tr_and,
                'log10_BF_TEP_flexible_vs_Null': l10_tf_null,
                'log10_BF_TEP_flexible_vs_Anderson': l10_tf_and,
                'bic_approximation_valid': {
                    'Anderson_vs_Null': valid_and_null,
                    'TEP_restricted_vs_Null': valid_tr_null,
                    'TEP_flexible_vs_Null': valid_tf_null,
                    'TEP_restricted_vs_Anderson': valid_tr_and,
                    'TEP_flexible_vs_Anderson': valid_tf_and,
                },
                'bic_approximation_note': (
                    "The BIC approximation BF ≈ exp(ΔBIC/2) is a large-sample result. "
                    "For n < 10 and extreme S/N it is unreliable.  When log10(BF) > 100 "
                    "the value is driven by formal uncertainties and should not be "
                    "reported as a literal probability ratio."
                ),
            },
            "model_selection": {
                'best_model_BIC': best_model_bic,
                'akaike_weights': akaike_weights,
                'interpretation': f"BIC selects {best_model_bic} model",
            },
        }

    def constant_floor_sensitivity_sweep(self, flybys, floors):
        """
        Transparency block: four-tier comparison under a range of constant
        shared systematic floors, showing the ordering dependence on sigma_sys
        (including the legacy anchor-dependent sigma_geom choice).  Each row
        records the constant floor and the resulting delta-BIC ordering.
        """
        rows = []
        for ss in floors:
            r = self._evaluate_four_tiers(flybys, float(ss), log_details=False)
            bics = r["information_criteria"]
            rows.append({
                "sigma_sys_const_mm_s": float(ss),
                "BIC_Null": bics["Null"]["BIC"],
                "BIC_Anderson": bics["Anderson"]["BIC"],
                "BIC_TEP_restricted": bics["TEP_restricted"]["BIC"],
                "BIC_TEP_flexible": bics["TEP_flexible"]["BIC"],
                "delta_BIC_TEP_restricted_vs_Null":
                    r["bayes_factors"]["delta_BIC_TEP_restricted_vs_Null"],
                "delta_BIC_TEP_restricted_vs_Anderson":
                    r["bayes_factors"]["delta_BIC_TEP_restricted_vs_Anderson"],
            })
        return rows

    def sign_agreement_model_comparison(self, eligible_flybys):
        """
        Sign-agreement-restricted four-tier comparison (NEAR, Galileo 1990, Rosetta 2005
        when Cassini is in the primary pool under relaxed sign gate).

        Emitted whenever at least two S/N-qualified flybys share the sign of the
        reference prediction at β_ref. When ``strict_sign_gate`` is true, this subset
        matches the primary gated ensemble.
        """
        sign_flybys = [
            fb
            for fb in eligible_flybys
            if flyby_sign_product(fb.get("dv_pred_base"), fb.get("dv_obs")) is not None
            and float(flyby_sign_product(fb["dv_pred_base"], fb["dv_obs"])) >= 0.0
        ]
        if len(sign_flybys) < 2:
            return None

        result = self.fractional_floor_model_comparison(sign_flybys)
        result["n_data"] = len(sign_flybys)
        result["included_flybys"] = [fb["name"] for fb in sign_flybys]
        return result


def main():
    """Execute stable model comparison."""
    logger = StepLogger("step_026_stable_model_comparison", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 026: FOUR-TIER MODEL COMPARISON")
    
    logger.section("FRAMEWORK")
    logger.info("Tiers: Null (0 param) | Anderson (2 param) | TEP restricted (1 param) | TEP flexible (3 param)")
    logger.info("TEP restricted is the primary evidence claim: all parameters except beta are pre-specified.")
    
    comparison = StableModelComparison(logger)
    
    # Load data
    loaded = comparison.load_data()
    if loaded is None:
        logger.error("Failed to load data")
        return 1

    flybys = loaded["eligible_flybys"]
    excluded = loaded["excluded_from_ensemble"]

    if not flybys:
        logger.error("No flybys passed Step 008 ensemble gates for model comparison")
        return 1

    logger.info(f"Loaded {len(flybys)} flybys passing Step 008 ensemble gates")
    logger.info(f"Catalog entries excluded from ensemble: {len(excluded)}")

    strict_sg = strict_sign_gate_from_config()
    logger.info(f"Ensemble config: strict_sign_gate={strict_sg} (pipeline default: false)")

    # Perform stable model comparison (primary gated ensemble)
    results = comparison.stable_model_comparison(flybys)
    results["strict_sign_gate"] = bool(strict_sg)
    results["ensemble_catalog"] = {
        "n_eligible": len(flybys),
        "excluded_from_ensemble": excluded,
        "ensemble_gate_policy": ENSEMBLE_GATE_POLICY,
    }

    sign_agreement = comparison.sign_agreement_model_comparison(flybys)
    if sign_agreement is not None:
        results["sign_agreement_model_comparison"] = sign_agreement
        logger.info(
            f"Sign-agreement diagnostic: n={sign_agreement['n_data']} "
            f"(TEP vs Null log10 B = "
            f"{sign_agreement['bayes_factors']['log10_BF_TEP_restricted_vs_Null']:.2f})"
        )

    # Ensemble composition robustness: test n=4 (Cassini included)
    n4_robustness = comparison.ensemble_composition_robustness(flybys, excluded, results)
    results["ensemble_composition_robustness"] = n4_robustness

    # Full-catalog breakthrough test: n=9 (all published flybys with explicit TEP predictions)
    loaded_full = comparison.load_full_catalog()
    if loaded_full and loaded_full["full_catalog_flybys"]:
        full_catalog = loaded_full["full_catalog_flybys"]
        full_results = comparison.full_catalog_model_comparison(full_catalog)
        results["full_catalog_model_comparison"] = full_results
        results["full_catalog_skipped"] = loaded_full["skipped"]
    else:
        logger.warning("Full-catalog comparison skipped: no flybys loaded")

    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step026_stable_model_comparison.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(convert_to_native_types(results), f, indent=2)
    
    logger.success(f"Stable model comparison complete. Saved to {output_file}")
    logger.add_output_file(output_file, "Stable model comparison results")
    
    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
