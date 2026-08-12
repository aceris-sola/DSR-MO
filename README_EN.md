# DSR-MO

> A single-objective niching algorithm for the **CEC 2026 Multimodal Optimization Benchmark** — **D**iverse-**S**earch & **R**efine for **M**ultimodal **O**ptimization
>
> Language: **[English](README_EN.md) | [简体中文](README.md)**

---

DSR-MO is a black-box, strategy-agnostic single-objective multimodal optimization algorithm. Its core idea is **"claim the landscape first, then cultivate it"**: the two hard sub-problems — *covering every peak basin* (discrete coverage) and *landing exactly on the basin floor* (continuous convergence) — are decoupled and handled by dedicated mechanisms, rather than squeezed simultaneously through a single objective signal.

> **🎯 Positioning: an extremely low-dimensional MMO algorithm**
>
> DSR-MO targets **D=2** specifically. The "claim first, cultivate later" philosophy relies on the **oracle-free self-estimated niche radius**, which is reliable only in low dimensions. In high dimensions (D≥5) the "curse of dimensionality" makes pairwise distances converge, degrading the self-estimated radius. DSR-MO therefore concentrates all its effort on the realistic, high-frequency D=2 regime rather than chasing unrealistic high-dimensional generalization.

---

## Table of Contents

- [1. Measured Results](#1-measured-results)
- [2. Problem We Address](#2-problem-we-address)
- [3. Algorithm Architecture](#3-algorithm-architecture)
- [4. Key Mechanisms & Math](#4-key-mechanisms--math)
- [5. Problems Encountered & Solutions](#5-problems-encountered--solutions)
- [6. Comparison with Others](#6-comparison-with-others)
- [7. Quick Start](#7-quick-start)
- [8. File Overview](#8-file-overview)
- [9. Reproduction](#9-reproduction)
- [10. License & Acknowledgements](#10-license--acknowledgements)

---

## 1. Measured Results

Protocol: **MATLAB R2026a, black box (official budget 20,000×D + oracle-free niche radius), DSR-MO final version (with reporting density filter)**, D=2 full 16 problems × 15 instances × 5 runs.

### 1.1 D=2 Summary

| Metric | Value |
|---|---|
| **Coverage (peak ratio / RPR)** | **0.8993** |
| **Precision** | **0.9921** |
| **F1** | **0.9345** |
| **SR (all-peaks success rate)** | 0.5383 |
| **Official score ½(RPR+F1)** | **0.9169** |

### 1.2 D=2 Per-Problem Results (mean over 15 instances × 5 runs)

| PID | Peaks | Coverage | Precision | F1 | PID | Peaks | Coverage | Precision | F1 |
|---|---|---|---|---|---|---|---|---|---|
| P01 | 20 | 0.966 | **1.000** | 0.982 | P09 | 10 | 0.960 | **1.000** | 0.979 |
| P02 | 20 | 0.987 | **1.000** | 0.993 | P10 | 10 | 0.997 | **1.000** | 0.999 |
| P03 | 20 | 0.971 | **1.000** | 0.985 | P11 | 10 | 0.993 | **1.000** | 0.996 |
| P04 | 20 | 0.954 | **1.000** | 0.976 | P12 | 10 | 0.979 | **1.000** | 0.989 |
| P05 | 20 | 0.995 | **1.000** | 0.998 | P13 | 10 | **1.000** | **1.000** | **1.000** |
| P06 | 20 | 0.617 | **1.000** | 0.758 | P14 | 10 | 0.532 | **1.000** | 0.684 |
| P07 | 20 | 0.872 | 0.942 | 0.904 | P15 | 10 | 0.612 | 0.932 | 0.734 |
| P08 | 20 | 0.975 | **1.000** | 0.987 | P16 | 10 | 0.977 | **1.000** | 0.988 |

- **Strong coverage**: P05/P10/P11/P12/P13 exceed 0.99, most reaching 1.000.
- **Precision ceiling**: after density filtering most problems reach Precision=1.000; only P07(~0.94)/P15(~0.93) are slightly lower.
- **Weakest**: P14(0.532)/P06(0.617)/P15(0.612) — the most deceptive, peak-dense problems where the self-estimated radius is limited; yet even there Precision ≈ 1.000.

## 2. Problem We Address

The core difficulty of multimodal optimization (MMO) is that **finding all peaks** and **finding them precisely** are hard to achieve simultaneously:

- **Finding all**: discover every equal-height global peak basin in the decision space (discrete coverage).
- **Finding precisely**: land exactly on the basin floor of each peak, without reporting redundant solutions (continuous convergence + compact reporting).

Classical approaches press both through a single objective signal, often sacrificing one for the other — strong convergence merges prematurely onto a few peaks, while strong spreading loses precision. DSR-MO's goal is to achieve **both high coverage and high precision without using true peak locations (oracle-free)**.

## 3. Algorithm Architecture

DSR-MO follows a three-stage **search → refine → report** pipeline. Each module is single-purpose and independently replaceable:

```
Initialize population
   │  (random sampling)
   ▼
Evolutionary search (main loop)
   │  variation (SBX crossover + polynomial mutation) → merge parents & offspring
   │  environmental selection: fitness + maximin farthest-point selection (keeps equal peaks)
   ▼
Peak extraction (one representative per basin)
   │  sort by objective + oracle-free niche-radius clearing
   ▼
Local refinement (derivative-free Nelder-Mead simplex)
   │  3× randomized restarts for D≥5
   ▼
Post-refinement objective threshold filter (discard local minima)
   ▼
Reporting-stage density filter (final report, oracle-free)
   │  simplified density clustering to deduplicate; keep the best per cluster
   ▼
Final solution set X, F
```

**Key design**: the search stage uses `maximin` to force the population to spread out in the decision space (preserving coverage), while the convergence accuracy is delegated entirely to the downstream Nelder-Mead refinement (preserving precision). The two are decoupled and do not interfere.

## 4. Key Mechanisms & Math

### 4.1 Environmental Selection: fitness + maximin

Each iteration merges parents and offspring, then selects `n` individuals back from the combined pool. The best-fit individual is kept first; then each selection picks the individual maximizing "better objective + farther from the already-selected set":

$$
\text{score}(x) \;=\; -\,\text{normF}(x)\;+\;\lambda\,\frac{d_{\min}(x)}{\mathrm{dom}}
$$

where:

- $\text{normF}(x) = \dfrac{F(x) - F_{\min}}{F_{\max} - F_{\min}}$ — normalized objective $\in[0,1]$, lower is better.
- $d_{\min}(x) = \min_{s\in S} \lVert x - s \rVert$ — distance from $x$ to the nearest already-selected point $S$.
- $\mathrm{dom} = \sqrt{\dim}\cdot(\mathrm{ub}-\mathrm{lb})$ — decision-space diagonal, used to normalize distances.
- $\lambda$ — diversity weight (default 20); larger $\lambda$ emphasizes decision-space footprint.

The two terms balance convergence ($-\text{normF}$) and diversity ($\lambda\, d_{\min}/\mathrm{dom}$), preserving all equal-height peaks.

### 4.2 Niche-based deduplication (clearing, oracle-free)

During peak extraction we sort by objective and classify a distant point as a new peak:

$$
d(p,\, \text{any existing representative}) > r \;\;\Longrightarrow\;\; p \text{ is a new peak representative}
$$

**Oracle-free radius** (does not rely on the true peak radius `problem.nichRad`):

$$
r \;=\; \operatorname*{median}_{i}\;\min_{j\neq i}\;\lVert x_i - x_j\rVert
$$

i.e. the median nearest-neighbor distance estimates the current landscape's peak spacing; adaptive and reproducible.

### 4.3 Local Refinement: Nelder-Mead derivative-free simplex

For each candidate peak representative, iterate Reflect → Expand → Contract → Shrink from an initial simplex, following curved valleys (e.g., Rosenbrock) to land exactly on a true global peak. For D≥5, 3 randomized restarts are performed per representative and the best result is kept, mitigating sensitivity to the initial simplex.

### 4.4 Reporting-stage density filter (final version)

A simplified density clustering (Union-Find connected components) is applied to the threshold-filtered report set; only the best representative (lowest objective) is kept per cluster, removing redundant reports:

$$
\varepsilon \;=\; \operatorname*{median}_{i}\;\min_{j\neq i}\;\lVert \hat{x}_i - \hat{x}_j\rVert,\qquad
\hat{x}_k \;=\; \operatorname*{arg\,min}_{x \in \text{cluster}_k} F(x)
$$

where $\varepsilon$ is the adaptive clustering radius (oracle-free); each cluster keeps the point with the lowest objective.

This module changes no search-core code and is an orthogonal, plug-and-play component — it is the key to Precision rising from 0.756 → 0.992.

### 4.5 CEC 2026 Evaluation Metrics

$$
\begin{aligned}
\text{RPR}\;=\;\frac{N_{\text{found}}}{N_{\text{GM}}} && \text{fraction of global peaks found}\\[2mm]
\text{Precision}\;=\;\text{RPR}\cdot\frac{N_{\text{GM}}}{N_{\text{sol}}} && \text{higher when reporting is compact}\\[2mm]
\text{F1}\;=\;\frac{2\,\text{Precision}\cdot\text{RPR}}{\text{Precision}+\text{RPR}}\\[2mm]
\text{Official}\ \text{Score}\;=\;\tfrac{1}{2}\big(\text{RPR}+\text{F1}\big)
\end{aligned}
$$

## 5. Problems Encountered & Solutions

During development and evaluation, DSR-MO encountered and solved the following problems one by one.

### 5.1 Porting bugs (benchmark correctness)

| Problem | Symptom | Solution |
|---|---|---|
| **Missing Weierstrass shift** | Missing `+0.5`; P07/P15 evaluated 1.8~5.2 off at true peaks | Added the shift; all 16 problems pass the benchmark self-check (maxDev=0) |
| **sequences.csv index out of range** | Official data is 1-based; an extra `+1` made `numUniform(10001)` out of range at D=5 | Fixed indexing; D=5 runs again |

### 5.2 Algorithm design problems

| Problem | Symptom | Solution |
|---|---|---|
| **Missing objective filter after refinement** | Ackley returned many local-minimum candidates, lowering Precision | Added a threshold filter; Precision improved markedly |
| **Over-budget + oracle radius** | Early runs used 50000×dim budget and true peak radius — not comparable to the official protocol | Switched to official budget 20000×dim + oracle-free radius |
| **Nelder-Mead high-dim sensitivity** | Sensitive to the initial simplex; unstable at D≥5 | 3 randomized restarts per representative at D≥5, keep best |
| **Redundant reports hurt precision** | The same basin reported repeatedly; Precision was only 0.756 | Added the reporting density filter (`dsr_mo_density_filter.m`); Precision rose to 0.992 with no loss of Coverage |

### 5.3 Methodological insight

Controlled experiments (isolating λ, budget, and oracle) reveal a counter-intuitive finding: **in the early stage of multimodal search, how far an individual is from the already-explored region matters far more than how good its objective is**. λ=20 forces the population to spread out; even if early objective values are poor, the downstream derivative-free Nelder-Mead refinement fully recovers the convergence loss. This is why "claim the landscape first, then cultivate it" works.

## 6. Comparison with Others

DSR-MO positions itself on the most mainstream track — the **IEEE CEC 2026 Multimodal Optimization Competition** — and openly challenges the year's **winner**. We do not shy away from this: you compete to go up against the best.

> Out of academic courtesy, we do not name the opponent here; "the CEC 2026 winner" refers to the winner algorithm of that competition (its paper: Chauhan & Dikshit, *S-CARD-CMSA: A Score-Aware Candidate Archive with Density-Filtered Reporting for Multimodal Optimization*, arXiv:2607.13764, 2026). All comparisons below are conducted in a respectful, fair, and same-protocol manner.

### 6.1 Measured Reproduction Comparison (D=2, main conclusion)

For fairness, DSR-MO reproduced the winner's open-source implementation on the **same machine (MATLAB R2026a), same benchmark**, running it problem-by-problem and instance-by-instance (16 problems × 15 instances, budget 20,000×D) for D=2, and scored every reported solution with the **exact same metric** used for DSR-MO (binary hit, relative tolerance accuracy 1e-3). Measured reproduction results (240 independent runs, one run per instance):

| Metric | Winner (reproduced in-house, D=2) | **DSR-MO** (D=2, 15 inst × 5 runs) |
|---|---|---|
| **Coverage / RPR** | **0.979** | 0.899 |
| **Precision** | 0.890 | **0.992** |
| **F1** | 0.916 | **0.935** |
| **Official score ½(RPR+F1)** | **0.947** | 0.917 |
| Success rate SR | 0.742 | 0.538 |

**Honest reading**: the in-house reproduction confirms the winner is highly competitive on D=2 — its **coverage (0.979)** clearly exceeds DSR-MO (0.899), a genuine strength of its search engine that DSR-MO freely acknowledges and respects. DSR-MO's edge lies in **report quality**: **Precision 0.992 and F1 0.935 both overtake** the winner, showing that DSR-MO's density-filtered deduplication delivers leaner, sharper final reports. The two form a complementary "strong-coverage search vs. precise-filter reporting" pattern — exactly the conclusion our controlled comparison is designed to surface.

> **Scope note**: the winner's numbers above are a one-run-per-instance in-house reproduction in the DSR-MO author's MATLAB environment (R2026a), scored with the same metric as DSR-MO; DSR-MO numbers are means over 15 instances × 5 runs. Both are D=2 with a 20,000×D budget and thus directly comparable. Reproduction script: `../复现_S-CARD-CMSA/run_reproduce_d2.m`.

**For reference**: the official numbers published in its technical report (RPR 0.805 / Precision 0.915 / F1 0.831 / official score 0.818) come from a 768-run **cross-dimension mixed** mean, not a pure D=2 full run — they are mentioned only for context and are **not** the basis for our conclusion. The in-house reproduction above is our comparison baseline.

### 6.2 Mechanism Comparison

| Dimension | Winner | DSR-MO |
|---|---|---|
| Search engine | Covariance-adaptive ES + taboo-region repulsion | Population evolution + fitness+maximin selection + Nelder-Mead refinement |
| Diversity | Active subpopulations + archive taboo repulsion of found peaks | maximin (farther from selected set is better) spreading |
| Convergence | Covariance-matrix-adaptive sampling | Derivative-free Nelder-Mead simplex |
| Reporting | Passive secondary archive + score-aware density filter | Threshold filter + simplified density clustering (oracle-free) |
| Uses true peaks | No | No |
| Core philosophy | Conservative: keep search, refine reporting | Active in both: strong coverage + report dedup |

**Conclusion**: based on the in-house reproduction (§6.1), DSR-MO overtakes the winner in the reporting (deduplication) stage — Precision 0.992 and F1 0.935 vs 0.890 / 0.916; DSR-MO trails slightly in the search (peak-finding) stage (coverage 0.899 < 0.979), which is honestly stated. Together they form a complementary "strong-coverage search vs. precise-filter reporting" pattern.

## 7. Quick Start

Requires **MATLAB R2026a** (or GNU Octave 11.3).

```matlab
addpath(pwd); addpath(fullfile(pwd, 'cec2026'));

% Full D=2: 16 problems × 15 instances × 5 runs, official budget 20000×dim
res = run_dsr_mo([2], 1:15, 1:16, 5, struct('pop_size', 300, 'lambda', 20));
```

Or run in one line from the command line:

```matlab
matlab -batch "addpath(pwd); addpath(fullfile(pwd,'cec2026')); res=run_dsr_mo([2],1:15,1:16,5,struct('pop_size',300,'lambda',20));"
```

## 8. File Overview

| File | Purpose |
|---|---|
| `dsr_mo.m` | Main algorithm (includes reporting density filter) |
| `dsr_mo_density_filter.m` | Reporting-stage density filter module |
| `dsr_mo_nelder_mead.m` | Local refinement (Nelder-Mead) |
| `dsr_mo_metrics.m` | Evaluation metrics (RPR / Precision / F1 / Coverage / SR) |
| `run_dsr_mo.m` | Main run entry |
| `make_figures.m` | One-click generation of all publication figures (output to `figs/`) |
| `figs/` | 6 evaluation figures (see §8.1) |
| `cec2026/` | CEC 2026 benchmark suite (`ProblemMM.m` / `BasicFun.m` / `UtilityMethod.m` / `pdist2.m` / `data/`) |
| `REPORT.md` | Detailed evaluation report |
| `README.md` / `README_EN.md` | This documentation (CN / EN) |

### 8.1 Figures (figs/)

| Figure | Content | Use |
|---|---|---|
| `fig1_bar_compare.png` | 16-problem Coverage & Precision grouped bars | Paper §4 / main figure |
| `fig2_radar.png` | DSR-MO vs CEC 2026 winner radar (4 public metrics) | Abstract / comparison |
| `fig3_scatter_peaks.png` | Decision-space peak distribution (P02/P14/P05) | Algorithm mechanism |
| `fig4_graft.png` | Density filter on/off comparison (P02/P05/P08) | §6 comparison evidence |
| `fig5_convergence.png` | Evaluation budget vs peaks found (P02) | Convergence analysis |
| `fig6_radius_box.png` | Self-estimated niche radius distribution (boxplot) | §4 parameter explanation |

## 9. Reproduction

```matlab
% MATLAB R2026a, official protocol: D=2 × 15 instances × 16 problems × 5 runs, black box (official budget + oracle-free)
% Final version: reporting density filter on by default (report_density=true), budget default 20000×dim
matlab -batch "addpath(pwd); addpath(fullfile(pwd,'cec2026')); res=run_dsr_mo([2],1:15,1:16,5,struct('pop_size',300,'lambda',20)); save(fullfile(pwd,'results_official.mat'),'res');"
```

## 10. License & Acknowledgements

- Benchmark: official CEC 2026 Multimodal Optimization Competition benchmark (`cec2026/`, by A. Ahrari et al.).
- Comparison target: the CEC 2026 Multimodal Optimization Competition winner (by Chauhan & Dikshit; open-source implementation in their repository).
- Algorithm code: DSR-MO (open-source, freely usable).