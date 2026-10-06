# Microscopic work under stationary boundary calibration

R. Arun Chandru | 6 October 2026.

## Extended abstract

A sufficiently high barrier is often described as an impenetrable wall. More accurate low-energy models retain an extrapolation length and further coefficients describing how the wall reflects a particle. These are stationary calibrations: they describe what the wall does when its microscopic structure is not being changed. This paper asks whether preserving such a calibration also fixes the work required to change the physical barrier.

The answer is no in the precise finite-potential models studied here. One can change a positive barrier along a nonconstant path while keeping any prescribed finite collection of its low-energy boundary-response coefficients exactly constant. The path and its amplitude are fixed before taking the high-wall limit. Its duration is proportional to the inverse wall height; the result does not describe a fixed laboratory duration as that height grows. Despite increasingly accurate agreement of stationary spectra and free energies, the drive produces rare energetic transitions. The paper proves a complete limiting law for changes of the ordered energy label, not merely a numerical observation of a few excited levels.

If the small wall scale is called `r`, label changes have probability proportional to `r^3`, and their characteristic energy is proportional to `r^-2`. Their mean dissipated work consequently decreases like `r`, while their work variance increases like `r^-1`. Most outcomes keep their energy label, but these typical outcomes need not have exactly zero work. All initially thermally populated modes and all final energy levels are included in the argument. Independent estimates control the probability, mean, and variance, so divergent moments are not inferred from an insufficient weak-convergence argument.

For an instantaneous profile change in a smooth bounded cavity, the leading mean and variance separate into a thermal boundary-pressure factor and microscopic profile quantities. This also holds for noninteracting Bose and Fermi gases with a fixed number of particles in the canonical ensemble. Statistics changes the absolute work coefficients. It cancels from one limiting normalized ratio; it does not make the full distributions identical, and the result is not about Bose condensation, interactions, or a thermodynamic limit.

An explicit smooth profile change provides a finite-wall test with exact scalar boundary calibration throughout its motion. Outward interval arithmetic encloses the physical finite-wall spectral and work quantities, including complete thermal tails. Positive dissipation and a bounded-energy detection event survive three specified positive durations under an explicit nondemolition reporting contract. A separate calculation compares six preparations, and an eight-profile design comparison improves a sufficient event-detection guarantee. These are deterministic theoretical certificates, not measurements from a laboratory experiment. The distinction between detecting an event and estimating its moments accurately is quantified: the latter can require far more preparations. The paper does not claim a resource-free sensing advantage or a demonstrated material implementation.


## Package contents

- `r_arun_chandru_calibrated_boundaries_microscopic_work.tex`: complete source, including bibliography and vector figure formulas/coordinates.
- `r_arun_chandru_calibrated_boundaries_microscopic_work.pdf`: compiled article, including appendices.
- `reproducibility/`: portable Python verification scripts and the three deterministic interval-certificate JSON files used for the reported numbers.
- `README.md`: this guide.


## Compile the article

Use an up-to-date TeX distribution containing REVTeX 4.2, AMS-LaTeX, `mathtools`, `booktabs`, `hyperref`, and PGFPlots. From the package directory:

```sh
pdflatex r_arun_chandru_calibrated_boundaries_microscopic_work.tex
pdflatex r_arun_chandru_calibrated_boundaries_microscopic_work.tex
```

Alternatively:

```sh
tectonic r_arun_chandru_calibrated_boundaries_microscopic_work.tex
```

Tectonic may download standard TeX packages on its first run. No shell escape or network access from manuscript code is required. The source uses an inline bibliography and inline vector plots, so neither BibTeX nor image generation is a separate compilation step.

## Reproduce the numerical certificates

Python scripts use only the standard library and same-directory imports. Python 3.12.14 was used for the verification runs. From `reproducibility/`, follow the sequential commands in its `README.md`. Run ordinary Python with `-B`, not optimized `-O`: assertions are part of the verification checks. Production enclosures use subdivision depth 12 (4096 cells). Smaller depth is a diagnostic only and must not be substituted for the reported numerical bounds.

Fresh runs of the smooth-wall, six-preparation, and eight-profile certificates, followed by the selected-candidate audit, reproduced every saved field exactly. Their measured runtimes on the verification machine were approximately 34, 38, 188, and 24 seconds, respectively; other machines may differ. These are interval enclosures with complete thermal-tail bounds, not Monte Carlo samples. The scripts and the mathematical derivations serve complementary roles; successful execution alone is not a proof of all analytical claims. The numerical intervals are outward enclosures, not confidence intervals. The sampling guarantees use separately specified measurement and independence assumptions.