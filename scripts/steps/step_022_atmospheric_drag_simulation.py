#!/usr/bin/env python3
"""
Step 022: Atmospheric Drag Simulation for Flyby Anomalies

Quantitatively tests whether atmospheric drag can explain observed flyby anomalies.
Computes atmospheric density at perigee altitudes using standard atmosphere models,
integrates drag force over hyperbolic trajectory, and compares to observed Δv.

This provides independent verification of the literature-cited exclusion
(atmospheric drag ~10^-6 mm/s at 1000-2000 km altitude).
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


@dataclass
class DragResult:
    """Results of atmospheric drag calculation for a single flyby."""
    mission: str
    perigee_altitude_km: float
    perigee_velocity_km_s: float
    atmospheric_density_kg_m3: float
    drag_acceleration_m_s2: float
    integrated_dv_mm_s: float
    observed_anomaly_mm_s: float
    drag_fraction_of_anomaly: float
    excluded: bool


class AtmosphericDragSimulator:
    """
    Simulates atmospheric drag effects on Earth flyby trajectories.

    Uses the standard piecewise-exponential atmosphere (US Standard
    Atmosphere 1976 base values, as tabulated in Vallado "Fundamentals of
    Astrodynamics and Applications" and Montenbruck & Gill "Satellite
    Orbits"): each layer i contributes rho(h) = rho_i * exp(-(h - h_i)/H_i)
    with layer-specific base density and scale height.  The thermospheric
    scale height grows with altitude (~50-270 km above 300 km); a single
    8.5 km scale height extrapolated above 100 km understates the density
    at 500-600 km by ~18 orders of magnitude and is only valid for the
    lower atmosphere.
    """

    # Layer table: (base altitude km, base density kg/m^3, scale height km)
    DENSITY_LAYERS = [
        (0.0,   1.225,     7.249),
        (25.0,  3.899e-2,  6.349),
        (30.0,  1.774e-2,  6.682),
        (40.0,  3.972e-3,  7.554),
        (50.0,  9.978e-4,  8.382),
        (60.0,  2.057e-4,  7.714),
        (70.0,  4.842e-5,  6.549),
        (80.0,  1.148e-5,  5.799),
        (90.0,  1.841e-6,  5.382),
        (100.0, 5.604e-7,  5.877),
        (110.0, 9.708e-8,  7.263),
        (120.0, 2.222e-8,  9.473),
        (130.0, 8.152e-9,  12.636),
        (140.0, 3.831e-9,  16.149),
        (150.0, 2.076e-9,  22.523),
        (180.0, 5.194e-10, 29.740),
        (200.0, 2.789e-10, 37.105),
        (250.0, 7.248e-11, 45.546),
        (300.0, 2.418e-11, 53.628),
        (350.0, 9.518e-12, 53.298),
        (400.0, 3.725e-12, 58.515),
        (450.0, 1.585e-12, 60.828),
        (500.0, 6.967e-13, 63.162),
        (600.0, 1.454e-13, 71.835),
        (700.0, 3.614e-14, 88.667),
        (800.0, 1.170e-14, 124.64),
        (900.0, 5.245e-15, 181.05),
        (1000.0, 3.019e-15, 268.00),
    ]

    def __init__(self):
        # Physical constants
        self.R_EARTH = 6378.137  # km (WGS84 equatorial radius)
        self.GM_EARTH = 3.986004418e14  # m^3 s^-2

        self.logger = StepLogger("step_022_atmospheric_drag_simulation")

    def atmospheric_density(self, altitude_km: float) -> float:
        """
        Compute atmospheric density at given altitude from the piecewise-
        exponential layer model.

        Parameters:
        -----------
        altitude_km : float
            Altitude above Earth's surface in km

        Returns:
        --------
        density_kg_m3 : float
            Atmospheric density in kg/m^3
        """
        altitude_km = max(0.0, altitude_km)
        # Select the highest layer whose base altitude <= altitude
        layer = self.DENSITY_LAYERS[0]
        for h0, rho0, H in self.DENSITY_LAYERS:
            if altitude_km >= h0:
                layer = (h0, rho0, H)
            else:
                break
        h0, rho0, H = layer
        return rho0 * np.exp(-(altitude_km - h0) / H)
    
    def drag_acceleration(
        self,
        density_kg_m3: float,
        velocity_m_s: float,
        area_m2: float = 10.0,
        drag_coefficient: float = 2.2
    ) -> float:
        """
        Compute drag acceleration: a_drag = 0.5 * ρ * v^2 * (C_d * A) / m
        
        Uses typical spacecraft parameters:
        - Area: 10 m^2 (representative cross-section)
        - Drag coefficient: 2.2 (typical for spacecraft)
        - Mass: 1000 kg (typical for interplanetary spacecraft)
        
        Parameters:
        -----------
        density_kg_m3 : float
            Atmospheric density
        velocity_m_s : float
            Spacecraft velocity
        area_m2 : float
            Cross-sectional area
        drag_coefficient : float
            Drag coefficient
        
        Returns:
        --------
        acceleration_m_s2 : float
            Drag acceleration in m/s^2
        """
        mass_kg = 1000.0  # Typical spacecraft mass
        return 0.5 * density_kg_m3 * velocity_m_s**2 * (drag_coefficient * area_m2) / mass_kg
    
    def integrate_drag_dv(
        self,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        drag_coefficient: float = 2.2,
        area_m2: float = 10.0,
        mass_kg: float = 1000.0,
    ) -> float:
        """
        Integrate drag acceleration along the osculating hyperbolic flyby
        trajectory to obtain the total velocity change.

        The perigee state (altitude, speed) fixes the osculating orbit:
        specific energy eps = v_p^2/2 - mu/r_p, angular momentum h = r_p v_p
        (the velocity is tangential at perigee), eccentricity
        e = sqrt(1 + 2 eps h^2/mu^2).  Parametrising by true anomaly nu,
        r(nu) = h^2/mu / (1 + e cos nu) and dt/dnu = r^2/h, so

            dv = int a_drag dt
               = int (1/2) rho(h(nu)) (C_d A/m) v(nu)^2 (r(nu)^2/h) dnu ,

        evaluated numerically over the full hyperbolic arc.  The integrand
        is exponentially concentrated around perigee, where the density is
        largest.

        Parameters:
        -----------
        perigee_altitude_km : float
            Perigee altitude
        perigee_velocity_km_s : float
            Perigee velocity
        drag_coefficient : float
            Spacecraft drag coefficient
        area_m2 : float
            Cross-sectional area
        mass_kg : float
            Spacecraft mass

        Returns:
        --------
        dv_mm_s : float
            Total velocity change from drag in mm/s
        """
        mu = self.GM_EARTH  # m^3 s^-2
        r_p = (self.R_EARTH + perigee_altitude_km) * 1000.0  # m
        v_p = perigee_velocity_km_s * 1000.0  # m/s
        K = (drag_coefficient * area_m2) / mass_kg  # m^2/kg

        eps = 0.5 * v_p * v_p - mu / r_p
        h_ang = r_p * v_p
        e = np.sqrt(1.0 + 2.0 * eps * h_ang * h_ang / (mu * mu))
        p = h_ang * h_ang / mu

        if e > 1.0:
            nu_lim = 0.999 * np.arccos(-1.0 / e)
        else:
            nu_lim = np.pi

        n_grid = 200001
        nu = np.linspace(-nu_lim, nu_lim, n_grid)
        r = p / (1.0 + e * np.cos(nu))
        alt_km = r / 1000.0 - self.R_EARTH
        rho = np.array([self.atmospheric_density(a) for a in alt_km])
        v2 = 2.0 * (eps + mu / r)
        # a_drag * dt/dnu = 0.5 * rho * K * v^2 * r^2 / h
        integrand = 0.5 * rho * K * v2 * r * r / h_ang
        dv_m_s = np.trapezoid(integrand, nu) if hasattr(np, "trapezoid") \
            else np.trapz(integrand, nu)
        return abs(dv_m_s) * 1000.0  # m/s -> mm/s
    
    def analyze_flyby(
        self,
        mission_name: str,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        observed_anomaly_mm_s: float
    ) -> DragResult:
        """
        Analyze atmospheric drag for a single flyby.
        
        Parameters:
        -----------
        mission_name : str
            Mission identifier
        perigee_altitude_km : float
            Perigee altitude
        perigee_velocity_km_s : float
            Perigee velocity
        observed_anomaly_mm_s : float
            Observed velocity anomaly (mm/s)
        
        Returns:
        --------
        result : DragResult
            Drag analysis results
        """
        # Compute atmospheric density
        rho = self.atmospheric_density(perigee_altitude_km)

        # Compute drag acceleration at perigee
        a_drag = self.drag_acceleration(rho, perigee_velocity_km_s * 1000)

        # Integrate along the osculating hyperbolic trajectory
        dv_drag = self.integrate_drag_dv(perigee_altitude_km, perigee_velocity_km_s)

        # Compare to observed anomaly where one is published. Atmospheric
        # drag removes kinetic energy, so its signed contribution to Delta_v_inf
        # is negative; the fraction is positive only when drag has the same
        # sign as the observed anomaly.
        if observed_anomaly_mm_s is not None and observed_anomaly_mm_s != 0:
            fraction = float(-dv_drag / observed_anomaly_mm_s)
            # Exclude if drag is < 1% of observed anomaly (wrong sign gives
            # a negative fraction and is likewise excluded)
            excluded = bool(fraction < 0.01)
        else:
            fraction = None
            excluded = None

        return DragResult(
            mission=mission_name,
            perigee_altitude_km=perigee_altitude_km,
            perigee_velocity_km_s=perigee_velocity_km_s,
            atmospheric_density_kg_m3=rho,
            drag_acceleration_m_s2=a_drag,
            integrated_dv_mm_s=dv_drag,
            observed_anomaly_mm_s=observed_anomaly_mm_s,
            drag_fraction_of_anomaly=fraction,
            excluded=excluded
        )
    
    def analyze_catalog(
        self,
        catalog_path: Path
    ) -> dict:
        """
        Analyze all flybys in the catalog.
        
        Parameters:
        -----------
        catalog_path : Path
            Path to step003_archival_flyby_catalog.json
        
        Returns:
        --------
        results : Dict
            Analysis results for all flybys
        """
        try:
            with open(catalog_path, 'r') as f:
                catalog = json.load(f)
        except (OSError, FileNotFoundError, json.JSONDecodeError) as e:
            self.logger.error(f"Failed to load catalog: {e}")
            return None
        
        results = []
        for flyby in catalog['flybys']:
            if not flyby['usable_for_analysis']:
                continue

            # Drag is evaluated for every usable flyby; the anomaly
            # fraction is only evaluated where a published anomaly exists.
            observed = flyby.get('published_anomaly_mm_s')
            if observed is not None and observed == 0:
                observed = None
            
            result = self.analyze_flyby(
                mission_name=flyby['mission_name'],
                perigee_altitude_km=flyby['perigee_altitude_km'],
                perigee_velocity_km_s=flyby['perigee_velocity_km_s'],
                observed_anomaly_mm_s=observed
            )
            results.append(result)
        
        # Summary statistics (fractions evaluated only on published anomalies)
        evaluated = [r for r in results if r.excluded is not None]
        all_excluded = all(r.excluded for r in evaluated)
        max_fraction = max((r.drag_fraction_of_anomaly for r in evaluated), default=None)
        dv_values = [r.integrated_dv_mm_s for r in results]

        return {
            'flyby_results': [
                {
                    'mission': r.mission,
                    'perigee_altitude_km': r.perigee_altitude_km,
                    'perigee_velocity_km_s': r.perigee_velocity_km_s,
                    'atmospheric_density_kg_m3': r.atmospheric_density_kg_m3,
                    'drag_acceleration_m_s2': r.drag_acceleration_m_s2,
                    'integrated_dv_mm_s': r.integrated_dv_mm_s,
                    'observed_anomaly_mm_s': r.observed_anomaly_mm_s,
                    'drag_fraction_of_anomaly': r.drag_fraction_of_anomaly,
                    'excluded': r.excluded
                }
                for r in results
            ],
            'summary': {
                'n_analyzed': len(results),
                'n_with_published_anomaly': len(evaluated),
                'all_excluded': all_excluded,
                'max_drag_fraction': max_fraction,
                'dv_range_mm_s': [min(dv_values), max(dv_values)] if dv_values else None,
                'conclusion': 'Atmospheric drag excluded for all flybys with published anomalies' if all_excluded else 'Atmospheric drag may contribute to some anomalies'
            }
        }


def main():
    """Execute atmospheric drag simulation."""
    logger = StepLogger("step_022_atmospheric_drag_simulation")
    
    try:
        simulator = AtmosphericDragSimulator()
        
        # Load flyby catalog
        catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
        
        if not catalog_path.exists():
            logger.error(f"Catalog not found: {catalog_path}")
            return None
        
        # Analyze all flybys
        results = simulator.analyze_catalog(catalog_path)
        
        # Save results
        output_path = PROJECT_ROOT / "results/step022_atmospheric_drag_simulation.json"
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Print summary
        print("\n" + "="*70)
        print("ATMOSPHERIC DRAG SIMULATION RESULTS")
        print("="*70)
        print(f"\nAnalyzed {results['summary']['n_analyzed']} flybys "
              f"({results['summary']['n_with_published_anomaly']} with published anomalies)")
        print("\nSummary:")
        print(f"  All excluded: {results['summary']['all_excluded']}")
        print(f"  Maximum drag fraction: {results['summary']['max_drag_fraction']:.2e}")
        print(f"  Δv range: {results['summary']['dv_range_mm_s']}")
        print(f"  Conclusion: {results['summary']['conclusion']}")

        print("\nPer-flyby results:")
        for r in results['flyby_results']:
            print(f"\n  {r['mission']}:")
            print(f"    Perigee altitude: {r['perigee_altitude_km']:.1f} km")
            print(f"    Atmospheric density: {r['atmospheric_density_kg_m3']:.2e} kg/m³")
            print(f"    Integrated Δv from drag: {r['integrated_dv_mm_s']:.2e} mm/s")
            if r['observed_anomaly_mm_s'] is not None:
                print(f"    Observed anomaly: {r['observed_anomaly_mm_s']:.2f} mm/s")
                print(f"    Drag fraction: {r['drag_fraction_of_anomaly']:.2e}")
                print(f"    Excluded: {r['excluded']}")
        
        if results['summary']['all_excluded']:
            logger.success(f"Atmospheric drag excluded for all {results['summary']['n_analyzed']} flybys")
        else:
            logger.info("Atmospheric drag not excluded for all flybys: " +
                        ", ".join(r['mission'] for r in results['flyby_results']
                                  if r['excluded'] is False) +
                        " (see summary conclusion)")
        logger.add_output_file(output_path, "Atmospheric drag simulation results")
        
        print(f"\n✓ Results saved to {output_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {e!s}")
        raise


if __name__ == "__main__":
    main()
