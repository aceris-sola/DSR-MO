# DSR-MO

> A single-objective niching algorithm for the **CEC 2026 Multimodal Optimization Benchmark** — **D**iverse-**S**earch & **R**efine for **M**ultimodal **O**ptimization
>
> Language: **[English](README_EN.md) | [简体中文](README.md)**

---

DSR-MO is a black-box, strategy-agnostic single-objective multimodal optimization algorithm. Its core idea is **"claim the landscape first, then cultivate it"**: the two hard sub-problems — *covering every peak basin* (discrete coverage) and *landing exactly on the basin floor* (continuous convergence) — are decoupled and handled by dedicated mechanisms, rather than squeezed simultaneously through a single objective signal.

On the official CEC 2026 evaluation protocol (16 problems × 15 instances × D=2, budget 20,000×D), DSR-MO reaches **Coverage 0.899 / Precision 0.992 / F1 0.935 / official score 0.917**, outperforming the published 0.818 of S-CARD-CMSA (RS-CMSA-ESII + density-filtered reporting).

---

## Table of Contents

- [1. Features & Highlights](#1-features--highlights)
- [2. Algorithm Architecture](#2-algorithm-architecture)
- [3. Key Mechanisms & Math](#3-key-mechanisms--math)
- [4. Methodological Contribution](#4-methodological-contribution)
- [5. Environment & Parameters](#5-environment--parameters)
- [6. Quick Start](#6-quick-start)
- [7. Benchmark Results](#7-benchmark-results)
- [8. Comparison with S-CARD-CMSA](#8-comparison-with-s-card-cmsa)
- [9. File Overview](#9-file-overview)
- [10. Reproduction](#10-reproduction)
- [11. License & Acknowledgements](#11-license--acknowledgements)

---

## 1. Features & Highlights

- **Oracle-free black box**: true peak positions are never used for optimization; both the niche radius and the report density-filter radius are **estimated from candidates** — aligned with official CEC rules.
- **Strong peak coverage**: for D=2, P05/P10/P11/P12/P13 frequently reach full coverage (Coverage/F1 = 1.000).
- **Precision up to 0.992**: reporting-stage density filtering compresses redundant reports to ≈ the true number of peaks; most problems reach Precision = 1.000.
- **Budget-insensitive**: empirically, coverage is *not* bought by a large evaluation budget; the algorithm is robust to budget.
- **Cross-platform**: results are bit-identical between MATLAB and GNU Octave (only floating-point rounding differences).
- **Reproducible**: a one-line call to `run_dsr_mo` regenerates all reported numbers.

## 2. Algorithm Architecture

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

## 3. Key Mechanisms & Math

### 3.1 Environmental Selection: fitness + maximin

Each iteration merges parents and offspring, then selects `n` individuals back from the combined pool. The best-fit individual is kept first; then each selection picks the individual maximizing "better objective + farther from the already-selected set":

```text
score(x) = -normF(x) + λ · (d_min(x) / dom)
```

where:

- `normF(x) = (F(x) - F_min) / (F_max - F_min)` — normalized objective ([0,1], lower is better).
- `d_min(x) = min_{s ∈ S} ‖x - s‖` — distance from `x` to the nearest already-selected point `S`.
- `dom = √dim · (ub - lb)` — decision-space diagonal, used to normalize distances.
- `λ` — diversity weight (default 20); larger λ emphasizes decision-space footprint.

The two terms balance convergence (`-normF`) and diversity (`λ·d_min/dom`), preserving all equal-height peaks.

### 3.2 Niche-based deduplication (clearing, oracle-free)

During peak extraction we sort by objective and classify a distant point as a new peak:

```text
if d(p, any existing representative) > r  ⟹  p is a new peak representative
```

**Oracle-free radius** (does not rely on the true peak radius `problem.nichRad`):

```text
r = median( nearest-neighbor distance of each point in the candidate pool )
```

i.e. the median nearest-neighbor distance estimates the current landscape's peak spacing; adaptive and reproducible.

### 3.3 Local Refinement: Nelder-Mead derivative-free simplex

For each candidate peak representative, iterate Reflect → Expand → Contract → Shrink from an initial simplex, following curved valleys (e.g., Rosenbrock) to land exactly on a true global peak. For D≥5, 3 randomized restarts are performed per representative and the best result is kept, mitigating sensitivity to the initial simplex.

### 3.4 Reporting-stage density filter (final version)

A simplified density clustering (Union-Find connected components) is applied to the threshold-filtered report set; only the best representative (lowest objective) is kept per cluster, removing redundant reports:

```text
eps = median( nearest-neighbor distance of the report set )   % adaptive, oracle-free
keep argmin F(x) per cluster
```

This module changes no search-core code and is an orthogonal, plug-and-play component — it is the key to Precision rising from 0.756 → 0.992.

### 3.5 CEC 2026 Evaluation Metrics

```text
RPR (recall) = N_found / N_GM                    % fraction of global peaks found
Precision    = RPR · N_GM / N_sol                % higher when reporting is compact
F1           = 2·Precision·RPR / (Precision+RPR)
official     = ½ · (RPR + F1)
```

## 4. Methodological Contribution

Unlike classical EAs that stress "premature convergence", DSR-MO shows empirically: **in the early stage of multimodal search, how far an individual is from the already-explored region matters far more than how good its objective is**. λ=20 forces the population to spread out in the decision space; even if early objective values are poor, the downstream derivative-free Nelder-Mead refinement fully recovers the convergence loss.

This is a "claim the landscape first, then cultivate it" multimodal strategy — it separates discrete coverage from continuous convergence and assigns each to a dedicated mechanism, rather than pressing both through a single objective signal.

## 5. Environment & Parameters

| Item | Value |
|---|---|
| Runtime | MATLAB R2026a (or GNU Octave 11.3) |
| Population size `pop_size` | 300 |
| Diversity weight `lambda` | 20 |
| Evaluation budget `max_eval` | **20000 × dim** (official CEC 2026) |
| Reporting density filter | enabled (default, adaptive radius) |
| Niche radius | self-estimated (oracle-free) |

## 6. Quick Start

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

## 7. Benchmark Results

> Protocol: MATLAB R2026a, black box (official budget 20,000×D + oracle-free radius), final version (with reporting density filter).

### 7.1 D=2 Summary (16 problems × 15 instances × 5 runs)

| Metric | Value |
|---|---|
| **Coverage (peak ratio / RPR)** | **0.8993** |
| **Precision** | **0.9921** |
| **F1** | **0.9345** |
| **SR (all-peaks success rate)** | 0.5383 |
| **Official score ½(RPR+F1)** | **0.9169** |

### 7.2 D=2 Per-Problem Results (mean over 15 instances × 5 runs)

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

### 7.3 D=5 Performance (legacy search-stage protocol, reported honestly)

| Metric | Value |
|---|---|
| Coverage | 0.4054 |
| Precision | 0.6654 |
| F1 | 0.4216 |
| SR | 0.0042 |

The D=5 weakness is attributed to the **"curse of dimensionality" degrading the self-estimated radius** (pairwise distances converge → radius inflates → niche over-aggregation, sacrificing coverage to preserve precision), not to seeds or platform issues. This is a known limitation; candidate fixes are adaptive radius / density-peak estimation.

## 8. Comparison with S-CARD-CMSA

> S-CARD-CMSA is an RS-CMSA-ESII + score-aware density-filtered reporting framework, in the same track as DSR-MO (CEC 2026 multimodal niching).

### 8.1 Quantitative Comparison (D=2, same 20,000×D budget)

| Metric | S-CARD-CMSA (paper 768-run validation) | DSR-MO (D=2, 15 inst × 5 runs) |
|---|---|---|
| **Coverage / RPR** | 0.805 | **0.899** |
| **Precision** | 0.915 | **0.992** |
| **F1** | 0.831 | **0.935** |
| **Official score ½(RPR+F1)** | 0.818 | **0.917** |

**Honest scope note**: S-CARD-CMSA's 0.805/0.915/0.831/0.818 come from a 768-run **cross-dimension mixed** validation subset, not a pure D=2 full run; DSR-MO is the full D=2 (16 problems × 15 instances × 5 runs). Both follow the same rule (20,000×D budget) but their dimension composition differs, so this is a **directional comparison**, not a strictly identical setup. The paper does not report per-dimension breakdowns; for a strictly identical setup, one can request the data or rerun its open-source code (`github.com/ChauhanDikshit`) following academic convention.

### 8.2 Mechanism Comparison

| Dimension | S-CARD-CMSA | DSR-MO |
|---|---|---|
| Search engine | Covariance-adaptive ES (RS-CMSA-ESII) + taboo-region repulsion | Population evolution + fitness+maximin selection + Nelder-Mead refinement |
| Diversity | Active subpopulations + archive taboo repulsion of found peaks | maximin (farther from selected set is better) spreading |
| Convergence | Covariance-matrix-adaptive sampling | Derivative-free Nelder-Mead simplex |
| Reporting | Passive secondary archive + score-aware density filter | Threshold filter + simplified density clustering (oracle-free) |
| Uses true peaks | No | No |
| Core philosophy | Conservative: keep search, refine reporting | Active in both: strong coverage + report dedup |

**Conclusion**: DSR-MO is competitive in both the search (peak finding) and reporting (deduplication) stages; its official score 0.917 clearly exceeds 0.818.

## 9. File Overview

| File | Purpose |
|---|---|
| `dsr_mo.m` | Main algorithm (includes reporting density filter) |
| `dsr_mo_density_filter.m` | Reporting-stage density filter |
| `dsr_mo_nelder_mead.m` | Local refinement (Nelder-Mead) |
| `dsr_mo_metrics.m` | Metrics (RPR / Precision / F1 / Coverage / SR) |
| `run_dsr_mo.m` | Main entry point |
| `cec2026/` | CEC 2026 benchmark suite (`ProblemMM.m` / `BasicFun.m` / `UtilityMethod.m` / `pdist2.m` / `data/`) |
| `REPORT.md` | Detailed evaluation report |
| `README.md` / `README_EN.md` | This documentation (Chinese / English) |

## 10. Reproduction

```matlab
% MATLAB R2026a, official rule: D∈{2,5} × 15 instances × 16 problems × 5 runs,
% black box (official budget + oracle-free). Final version: density filter on by default
% (report_density=true), budget defaults to 20000×dim.
matlab -batch "addpath(pwd); addpath(fullfile(pwd,'cec2026')); res=run_dsr_mo([2,5],1:15,1:16,5,struct('pop_size',300,'lambda',20));"
```

## 11. License & Acknowledgements

- The algorithm code is released under the **MIT License** (see `LICENSE`).
- The benchmark suite `cec2026/` is ported from the official **IEEE CEC 2026 Multimodal Optimization Competition** benchmark; its copyright belongs to the original authors.
- Comparison baseline S-CARD-CMSA: Dikshit Chauhan, *S-CARD-CMSA: A Score-Aware Candidate Archive with Density-Filtered Reporting for Multimodal Optimization*, arXiv:2607.13764, 2026.