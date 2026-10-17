# The Temporal Equivalence Principle and the Earth Flyby Anomaly: Clock-Sector Response and Coherent-Link Cancellation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454862.svg)](https://doi.org/10.5281/zenodo.19454862)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
![TEP-EFA: Earth Flyby Anomaly](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.4 (Yogyakarta)  
**First published:** 17 May 2026 · **Last updated:** 17 October 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862)  
**Website:** [https://mlsmawfield.com/tep/efa/](https://mlsmawfield.com/tep/efa/)  
**Paper Series:** TEP Series · Group C (III. Laser Ranging, Spacecraft Tracking & Solar-System Tests) · Paper 15 (`C15`)

## Abstract

<!-- Generated from site/components/1_abstract.html; update the manuscript source, then run tools/sync_current_corpus_readmes.py. -->

Twelve Earth gravity assist flybys spanning nine spacecraft are analyzed within the Temporal Equivalence Principle (TEP) framework. TEP posits that proper time is governed by a dynamical scalar field φ; all non-gravitational matter couples universally to a causal matter metric through the conformal factor A(φ) = exp(β_A φ/M_Pl) with β_A = −1. Under Paper 0's amplitude–shear split (v0.16), the Earth-vicinity Temporal Shear is strongly screened by the two-body kinetic operator, so the flyby anomaly cannot be a shear-force effect. The candidate mechanism is instead the clock sector: highly asymmetric, fast hyperbolic flybys accumulate rapid proper-time offsets at perigee through the unsuppressed amplitude channel S_A, which orbit-determination codes — reconstructing velocities under an assumed GR time standard — misread as ΔV discontinuities wherever the tracking link references the spacecraft's own time standard. The channel analysis performed here shows the term survives only in clock-carrying observables (one-way, non-coherent, or regenerative-ranging links and onboard-clock telemetry); in the coherent-transponder class assigned to the published catalogue from the literature, the transponder endpoint cancels it identically, so on that classification the catalogue anomalies are not produced by the clock sector in the measured observable (encounter-specific tracking modes remain to be checked) — the mechanism's claim is instead a falsifiable clock-carrying prediction whose m/s-scale raw amplitude is conditional on the adopted profile and the eccentric-orbit clock constraint in Section 5.7a. The fitted response profile S_⊕(r) is accordingly reinterpreted as the clock-sector amplitude calibration, with the characteristic terrestrial field scale R_sol = R_T(M_⊕) = (3M_⊕/4πρ_T)¹/3 ≈ 4146 km for ρ_T ≈ 20 g/cm³. Since the UCD calibration of ρ_T descends from the GNSS coherence-length measurement, this terrestrial length scale is a single shared input rather than an independent validation.

The accumulated proper-time offset manifests as an apparent velocity discontinuity in GR-based orbit reconstruction — a clock-sector artefact misread as ΔV without an unmodeled force — in the measurement classes that carry the spacecraft's own time standard. The non-radial component of the offset integral, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, the velocity-dependent response template, and the perigee plasma template layer, supplies the geometric template against which the observed flyby anomaly pattern is organized. Six published nonzero entries are retained in the re-transcribed catalogue — four firm detections (NEAR, Galileo 1990, Rosetta 2005, Galileo 1992) and two marginal ~2σ entries (Cassini, MESSENGER) — alongside three published nulls/bounds and three flybys without public anomaly reports; the re-transcription against the Anderson et al. (2008) primary source corrected several sign, uncertainty, and declination entries (Section 2.3). Under a full-catalog likelihood with a fractional misspecification floor — a per-flyby systematic set to the measured between-flyby amplitude heterogeneity ($CV_\beta \approx 0.47$) times each model's own predicted signal, iterated to a fixed point and invariant under the arbitrary field-normalization anchor — the archived information-criterion comparison describes how the restricted geometric template organizes the catalogue relative to the Null by $\Delta{\rm BIC} \approx 1.67\times10^{4}$; the Anderson empirical baseline, whose two parameters are fitted on the same catalogue, attains a nominally smaller BIC but at parity-level separation ($\Delta{\rm BIC} \approx 4.8$), an ordering that reverses in TEP restricted's favour ($\Delta{\rm BIC} \approx +2767$) when the published uncertainties are taken at face value without a systematic floor. These scores describe empirical geometry-template organization, not evidence for the clock-sector mechanism, which is absent from the coherent catalogue observable. They retain the archived iteratively weighted fitting convention; a fully optimized prediction-dependent-variance likelihood and its nuisance accounting are required before interpreting them as calibrated model-selection evidence.

A structural consistency check is provided by the ordering of anomaly magnitudes against the trajectory-asymmetry factor (cosδ_in − cosδ_out). In the gated n = 6 detection ensemble the ordering is positive but not decisive (Spearman ρ = +0.66, p = 0.156): NEAR carries both the largest asymmetry factor and the largest published anomaly in every ranked subset. Across the broader catalogue the rank association is moderate: Spearman ρ = 0.56 (p = 0.11) over the n = 9 flybys with published asymmetries, and ρ = 0.68 (p = 0.090) over the n = 7 flybys with nonzero published anomalies, computed in the pipeline's ordering diagnostic. These are reported as a qualitative consistency check rather than a primary discriminator. The inverse-variance weighted mean $\beta_{\rm fit} \approx 4.74 \times 10^{-3}$ (formal SE $5.83 \times 10^{-4}$; random-effects SE $1.28 \times 10^{-3}$; between-flyby $\tau \approx 2.23 \times 10^{-3}$) quantifies the cross-flyby effective response scale; it is not the bare conformal coupling $\beta_A$. The amplitude-informative per-flyby fits span $4.43 \times 10^{-3}$ to $1.44 \times 10^{-2}$ (with the negative-branch Galileo 1992 member at $6.7 \times 10^{-2}$), consistent with geometry-dependent modulation. The fitted amplitudes are trajectory-response coefficients, not the Cassini source charge; Cassini PPN compliance is evaluated through the solar source-charge projection $\alpha_{\rm eff} = S_\Sigma \alpha_0$ (Section 4.6.1a). The re-transcribed catalogue resolves the historical Cassini sign question: the primary-source anomaly is negative (−2 ± 1 mm/s), so the kernel's predicted sign agrees and the residual tension is magnitude-only. The geometry envelope carries seven nominal deterministic coefficients plus one fitted amplitude against six retained nonzero entries, so the model-selection comparisons are read as organisational support under the adopted convention rather than calibrated evidence.

This work shows that a restricted trajectory-geometry template quantitatively organizes the published Earth flyby anomaly catalogue under the adopted Anderson trajectory-geometry convention, while remaining consistent with precision solar system constraints and identifying specific stress tests—principally the Galileo 1992 amplitude shortfall and the Juno raw-tension null, together with raw DSN reanalysis—that are required to advance from model organisation to key confirmation. The catalogue's evidential status is conditional on one provenance question carried openly here: the orbit-determination filtering interpretation invoked for the published nulls is symmetric with the alternative that the early anomalies were artefacts of less sophisticated orbit-determination pipelines — the two readings are equally consistent with the present data, and only an independent minimal-OD reanalysis of anomaly and null flybys under identical force models can arbitrate between them.

The terrestrial response profile $S_\oplus(r)$ constitutes the Earth-vicinity geometric realization of the clock-sector amplitude response $S_A(\mathcal{E})$. Rather than a hard boundary, this radial profile traces the continuous amplitude structure of the temporal field within Earth's local potential well.

Keywords: Earth flyby anomaly, Temporal Equivalence Principle, clock-sector artefact, proper-time offset, trajectory asymmetry, geometric screening, Temporal Topology, Temporal Shear

Current manuscript source: [1_abstract.html](site/components/1_abstract.html).

## Key Findings

The current findings, numerical results and comparison assumptions are reported in [Results](site/components/5_results.html), [Conclusions](site/components/7_conclusion.html). The abstract above is synchronized from the current manuscript.

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 15** | **TEP-EFA** (This repo) | Temporal Equivalence Principle and the Earth Flyby Anomaly: Clock-Sector Response and Coherent-Link Cancellation | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454619](https://doi.org/10.5281/zenodo.19454619) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |

## Directory Structure

```
TEP-EFA/
├── data/                          # Flyby trajectories and processed outputs
│   └── raw/                       # Original data sources
├── scripts/
│   ├── steps/                     # Sequential analysis pipeline
│   └── utils/                     # Shared utilities
├── results/                       # Analytical outputs and figures
├── site/
│   ├── components/                # Manuscript HTML sections
│   └── public/                    # Static assets
├── logs/                          # Pipeline execution logs
├── config/                        # Pipeline configuration
├── README.md                      # This file
└── requirements.txt               # Python dependencies
```

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-EFA.git
cd TEP-EFA

# Install dependencies
pip install -r requirements.txt
```

## Reproduction Steps

The flyby TEP analysis pipeline is fully automated. Run the complete pipeline:

```bash
python scripts/run_all.py
```

This will populate `results/` and `data/processed/` with fresh analysis outputs.

## Citation

```bibtex
@article{smawfield2026earthflyby,
  title={The Temporal Equivalence Principle and the Earth Flyby Anomaly: Clock-Sector Response and Coherent-Link Cancellation},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.19454862},
  note={Preprint v0.4 (Yogyakarta)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
