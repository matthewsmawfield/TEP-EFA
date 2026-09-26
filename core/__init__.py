"""TEP Core — canonical Python package for the Temporal Equivalence Principle.

Version: TEP v0.14 (Jakarta)

This package provides the shared physics layer used by all TEP papers:
  - constants: physical and phenomenological parameters
  - conformal_scaling: A(phi), Temporal Shear, effective coupling
  - cosmology: LCDM + TEP distance-redshift relations
  - scalar_field: lab-scale scalar-field solver
  - screening: environment-dependent screening functions
  - evidence: Bayesian evidence and model-comparison utilities
"""

from . import conformal_scaling, constants, cosmology, evidence, scalar_field, screening

__all__ = [
    "conformal_scaling",
    "constants",
    "cosmology",
    "evidence",
    "scalar_field",
    "screening",
]

# Re-export commonly used symbols at package level for convenience
from .conformal_scaling import (
    conformal_factor,
    conformal_factor_small,
    effective_g,
    g_eff_variance,
    minimum_steepness_for_retention,
    minimum_steepness_for_suppression,
    screening_diagnostics,
    temporal_shear_from_scalar_field,
)
from .constants import (
    ALPHA_LOG,
    BETA_A,
    BETA_GEOM,
    C_LIGHT,
    G_CM3_TO_KG_M3,
    G_NEWTON,
    GNSS_LAMBDA_T_EXPONENTIAL_BY_CENTER,
    GNSS_LAMBDA_T_LONGSPAN_CODE_ERR_KM,
    GNSS_LAMBDA_T_LONGSPAN_CODE_KM,
    ILLUSTRATIVE_BETA_A,
    KG_M3_TO_G_CM3,
    LAB_COHERENCE_LENGTH_M,
    LAMBDA_T_MGEX_ERR_KM,
    LAMBDA_T_MGEX_KM,
    LAMBDA_T_MGEX_R2,
    M_PLANCK,
    M_REF,
    M_SUN,
    MPC_TO_M,
    RHO_C,
    RHO_T,
    SCREENING_LENGTH_KM,
    VERSION,
    VERSION_CODENAME,
    VERSION_STRING,
)
from .scalar_field import (
    compute_temporal_shear_from_mass_gradient,
    scalar_field_difference,
    scalar_field_logarithmic,
    solve_scalar_field_cylinder,
    solve_scalar_field_layered,
    solve_scalar_field_layered_weighted,
    solve_scalar_field_uniform_density,
)
from .screening import (
    coupling_screening_factor,
    matter_acceleration,
    physical_shear,
    screening_factor,
    screening_ratio,
    universal_screening_function,
)
