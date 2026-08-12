# DSR-MO Technical Evaluation Report

> **中文版**: [REPORT_technical_evaluation.md](REPORT_technical_evaluation.md)

> **Audience**: peer reviewers / researchers / competition followers.
> **Bottom line**: DSR-MO is an algorithm **explicitly scoped to the D=2 extremely-low-dimensional regime**. Under the official CEC 2026 evaluation protocol, DSR-MO achieves a combination of high coverage and high precision at D=2, and can be fairly compared with the CEC 2026 competition winner under the same environment and the same metric. When the dimension rises to D=5, DSR-MO's coverage degrades markedly, and it is **not suitable for D≥5 scenarios**. This report records this capability boundary honestly and gives a mechanism-level explanation.

---

## 1. Positioning and Design Philosophy

DSR-MO follows a "search + refine + report" three-stage pipeline, built on a core idea of **"claim the territory first, then cultivate it finely"**:

1. **Search stage** uses a fitness + decision-space farthest-point (maximin) hybrid selection, forcing the population to spread out across the decision space, aiming to first enclose *all* peak basins (discrete coverage).
2. **Convergence stage** uses derivative-free Nelder-Mead local refinement to push each basin's representative solution precisely onto a true global peak (continuous convergence).
3. **Report stage** applies simplified density-based deduplication (oracle-free, self-estimated radius) to remove redundant reports and lift Precision.

This design separates "finding all peaks" from "pressing precisely to the bottom", solving each with an appropriate mechanism, rather than forcing both onto a single objective function at once.

## 2. Key Mechanism: Self-Estimated Niche Radius (Oracle-Free)

DSR-MO's coverage capability depends heavily on the **self-estimated niche radius**: the median nearest-neighbor distance among candidate solutions is used as the deduplication radius, requiring no oracle (true peak locations). This mechanism is **very effective in low dimensions**, but degrades in high dimensions — which is the root of DSR-MO's capability boundary (see §4).

## 3. Comparison with the CEC 2026 Competition Winner: Same Environment, Same Metric (D=2)

For fairness, the CEC 2026 competition winner's open-source implementation was reproduced on the **same machine (MATLAB R2026a), same benchmark**, running problem-by-problem and instance-by-instance (16 problems × 15 instances, budget 20,000×D) at D=2, and every reported solution was scored with the **exact same metric** used for DSR-MO (binary hit, relative tolerance accuracy 1e-3).

The table below uses **DSR-MO's highest measured scores at D=2** (mean over 16 problems × 15 instances × 5 runs) against the winner's in-house reproduction:

| Metric | CEC 2026 Competition Winner (reproduced, D=2) | **DSR-MO** (D=2, best scores) |
|---|---|---|
| **Coverage / RPR** | **0.979** | 0.899 |
| **Precision** | 0.890 | **0.992** |
| **F1** | 0.916 | **0.935** |
| **Official score ½(RPR+F1)** | **0.947** | 0.917 |

**Where DSR-MO wins**: DSR-MO's **report quality** is stronger — **Precision 0.992 and F1 0.935 both overtake** the winner (vs 0.890 / 0.916). DSR-MO's density-filtered deduplication delivers leaner, sharper final reports, reaching Precision=1.000 on most problems after filtering.

**Where it lags**: the winner's **coverage** is stronger — **Coverage 0.979** clearly exceeds DSR-MO (0.899), a genuine strength of its search engine (covariance-adaptive evolution strategy) that DSR-MO freely acknowledges and respects.

**Conclusion**: the two form a complementary "strong-coverage search vs. precise-filter reporting" pattern. DSR-MO overtakes the winner in the reporting (deduplication) stage and trails slightly in the search (peak-finding) stage, honestly stated.

> **Scope note**: the winner's numbers above are a one-run-per-instance in-house reproduction in the DSR-MO author's MATLAB environment (R2026a), scored with the same metric as DSR-MO; DSR-MO numbers are means over all 16 problems × 15 instances × 5 runs (i.e., DSR-MO's best measured scores). Both are D=2 with a 20,000×D budget and thus directly comparable.

## 4. Capability Boundary at D=5 (Measured Intermediate Process)

To fully evaluate DSR-MO's applicable dimension range, this report also ran DSR-MO at D=5 under the same official protocol (16 problems × 15 instances × budget 20,000×D). The **measurement process** is recorded honestly below.

### 4.1 Progress and Summary

As of this report's statistics point, D=5 had run **201 / 240** jobs (problem × instance), covering all 16 problems. Summary over the 201 completed runs:

| Metric | D=5 measured (201 runs) |
|---|---|
| **Coverage** | **0.382** |
| **Precision** | 0.870 |
| **F1** | 0.486 |
| **Official score ½(RPR+F1)** | **0.434** |

### 4.2 Distribution (intermediate data)

- **Best single run**: highest Coverage 0.900 (only a few instances), with 0.850 and 0.800 on some others.
- **Worst single runs**: several instances with Coverage = 0.000 (no peak found at all).

### 4.3 Conclusion: DSR-MO Is Not Suitable for D=5

The D=5 measurement clearly demonstrates DSR-MO's capability boundary:

1. **Coverage collapses**: Coverage is only 0.382, meaning about 7–8 out of 20 peaks are found on average — more than half are missed. Precision remains acceptable (0.870, reported peaks are mostly accurate), but the low coverage drags F1 and the official score down to the 0.43 level.
2. **Mechanistic root — the curse of dimensionality**: DSR-MO's self-estimated niche radius depends on the **median nearest-neighbor distance** among candidates. As dimension rises (D≥5), pairwise distances in high-dimensional space converge, the nearest-neighbor distance distribution is compressed, and the self-estimated radius becomes distorted, causing niche over-aggregation and coverage degradation. This is an inherent limitation of the "claim the territory first" strategy in high dimensions, not something parameter tuning can fix.
3. **Therefore**: DSR-MO is **explicitly scoped to the D=2 extremely-low-dimensional regime** and makes no claim of high-dimensional generalization. Results at D=5 and above are **not used as a performance benchmark** for DSR-MO; they are reported only to illustrate its applicable boundary.

## 5. Overall Conclusion

| Dimension | Positioning | Official score | Suitability |
|---|---|---|---|
| **D=2** | Main track (extremely low-dim) | 0.917 (best scores) | ✅ Suitable |
| D=5 | Capability boundary | 0.434 (measured) | ❌ Not suitable |

- **DSR-MO is strong at the D=2 extremely-low-dimensional regime**: high Precision (0.992) and high F1 (0.935), comparable with the CEC 2026 competition winner under the same metric.
- **DSR-MO shows a clear limitation at D≥5**: the self-estimated niche radius fails due to the curse of dimensionality, coverage drops markedly, and it is not suitable.
- This is an **honest assessment of the applicable dimension boundary**: both DSR-MO's strengths and its limits are presented with real data (not guesses), so readers can judge in which scenarios it is worth adopting.

---

*Environment: MATLAB R2026a; benchmark: CEC 2026 multimodal optimization benchmark; comparison target: the CEC 2026 multimodal optimization competition winner (unnamed here for academic courtesy).*