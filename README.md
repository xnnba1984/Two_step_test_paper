This repository contains the analysis and simulation code that produces all empirical results in the manuscript "_Xi N.M., Wang L. (2025). A Two-Step Test for Zero-Inflated Biomarkers in Early-Phase Clinical Trials_." Each script is independent of the ```twostepAKSA``` R package and no external functions are required. The preprint can be downloaded [here]().

## Repository overview
### aksa_power_loss.R
Demonstrates raw AKSA loses power under tail‑only alternatives as the zero‑inflation rate increases. The script simulates outcomes under a monotone treatment effect confined to the positive biomarker tail, evaluates AKSA power across grids of sample size and zero-finlation rate π₀, and produces the power‑loss panel used in the manuscript.

### grid_null.R
Type‑I error study under the global null with independent component tests. Builds a factorial grid over sample size and zero‑inflation, simulates datasets with no treatment effect, computes p‑values for the spike test (A), tail test (B), raw AKSA, and their Fisher/Brown combinations, and summarizes empirical level at α = 0.05.

### grid_spike.R
Power study for the spike‑only alternative (constant treatment benefit among subjects with X = 0; no benefit in the positive tail). Runs the same factorial grid as above, estimates power for AKSA and for the two‑step test with Fisher and Brown combinations, and saves summaries/figures.

### grid_mono.R
Power study for the tail‑only alternative (linearly increasing benefit with biomarker rank among X > 0; no effect at the spike). Mirrors grid_spike.R but with monotone effects.

### grid_mix.R
Power study for the mixed alternative (constant spike effect + monotone tail effect). Uses the same grid and output structure as the other grid_* scripts.

### cor_null.R
Type‑I error under positively correlated component tests. Generates pairs of correlated Uniform(0,1) p‑values by transforming bivariate Normals with correlation ρ and evaluates Fisher vs Brown across ρ ∈ {0, 0.1, …, 0.8}. The script also records a correlation‑specific Fisher cutoff (critical p‑value) that attains 0.05 level under the dependent null that is used downstream by the power simulations.

### cor_spike.R
Power under spike‑only alternatives when the spike and tail p‑values are positively correlated. Induces controlled dependence between spike and tail test by adding a common random perturbation to the treated outcomes, applies Brown at nominal α = 0.05 and Fisher at the ρ‑adjusted cutoff from cor_null.R, and reports Brown–Fisher power comparisons. 

### cor_mono.R
Same as cor_spike.R but for the tail‑only alternative. 

### cor_mix.R
Same as cor_spike.R but for the mixed alternative. 

### skew_null.R
Type‑I error under skewed positive‑tail biomarker distributions. Replaces the Uniform(0,1) tail by right‑skewed Beta distributions (mild, moderate, strong) while preserving zero‑inflation and randomization; evaluates empirical type-I errors for AKSA, Fisher, and Brown across π₀ and N. 

### skew_spike.R
Power under spike‑only alternatives with right‑skewed positive tails. Same factorial grid and outputs as the independent‑tail spike setting. 

### skew_mono.R
Power under tail‑only alternatives with skewed positive tails. Produces the monotone/skew panels in the manuscript. 

### skew_mix.R
Power under mixed alternatives with skewed positive tails. Complements the spike and monotone skew analyses. 

## Prerequisites
R ≥ 4.1

Minimal dependencies across scripts are: Rcpp, RcppArmadillo (for the fast two‑sample KS), ggplot2, dplyr, tidyr, viridisLite

Install required R package with
```r
install.packages(c("Rcpp","RcppArmadillo","ggplot2","dplyr","tidyr","viridisLite"))
```

## Run the analyses

Independent tests: ```grid_null.R```, ```grid_spike.R```, ```grid_mono.R```, ```grid_mix.R```.

Dependence calibration: run ```cor_null.R``` once to produce the Fisher correlation‑specific cutoff.

Power under test correlation: ```cor_spike.R```, ```cor_mono.R```, ```cor_mix.R```.

Right‑skew sensitivity: ```skew_null.R```, ```skew_spike.R```, ```skew_mono.R```, ```skew_mix.R```.

Illustration of AKSA’s tail‑only loss: ```aksa_power_loss.R```.

Each script constructs its own factorial grid (sample size, π₀, effect sizes), sets a reproducible seed, runs Monte Carlo replicates, and writes summaries/figures to disk.

## Computation
The default settings (1,000 permutations × 1,000 replicates per grid point) can be time‑consuming. For a quick test, temporarily reduce the number of permutations/replicates near the top of each script, then restore defaults for final figures.


## Contact
For any questions or issue reports, pleasae open an issue or email nxi@ucla.edu.
