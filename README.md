# DSR-MO

**DSR-MO**（Diverse-Search & Refine for Multimodal Optimization）是一个用于 **CEC 2026 多模态优化基准** 的单目标 niching 算法。

核心思想是 **"先占地盘、后精耕细作"**：把"找全所有峰盆地"（离散覆盖）与"精确压到盆底"（连续收敛）两个难度分离，分别交给不同机制，而不是在单一目标函数上同时硬压两者。

## 算法机制

1. **fitness + 最远点（maximin）混合环境选择**：先保留目标最优个体，之后每次选「目标越好 + 离已选集合越远」得分最高的个体，在收敛性与决策空间足迹之间取得平衡，从而保留全部等高峰。
2. **生态位去重（clearing）**：按**自估半径**（候选解最近邻距离中位数，无先知信息）每盆地提取一个代表解。
3. **Nelder-Mead 局部精修**：对候选峰做无导数单纯形搜索，沿弯曲谷道精确压到真全局峰；D≥5 时对每个代表做 3 次随机扰动重启，缓解对初始单纯形的敏感性。
4. **精修后函数值阈值过滤**：剔除局部极小候选，提升精度。
5. **报告期密度过滤**（`dsr_mo_density_filter.m`）：对报告集做 DBSCAN 式密度聚类去冗余——同一盆地内的重复代表只保留最优者，消除"多报冗余"导致的 Precision 损失。半径自适应（报告集最近邻距离中位数），全程无先知。

## 文件清单

| 文件 | 作用 |
|---|---|
| `dsr_mo.m` | 主算法（含报告期密度过滤） |
| `dsr_mo_density_filter.m` | 报告期密度过滤模块 |
| `dsr_mo_nelder_mead.m` | 局部精修 |
| `dsr_mo_metrics.m` | 评测指标（RPR / Precision / F1 / Coverage / SR） |
| `run_dsr_mo.m` | 主运行入口 |
| `cec2026/` | CEC2026 基准套件（移植自官方 Python 版） |
| `REPORT.md` | 详细评测报告 |

## 运行

需要 **MATLAB R2026a**（或 GNU Octave 11.3，见 `cec2026/` 兼容说明）。

```matlab
addpath(pwd); addpath(fullfile(pwd, 'cec2026'));
% D=2 全量：16 题 × 15 实例 × 5 次，官方预算 20000×dim
res = run_dsr_mo([2], 1:15, 1:16, 5, struct('pop_size', 300, 'lambda', 20));
```

或命令行一键：

```matlab
matlab -batch "addpath(pwd); addpath(fullfile(pwd,'cec2026')); res=run_dsr_mo([2],1:15,1:16,5,struct('pop_size',300,'lambda',20));"
```

## 评测结果（D=2 全量，MATLAB R2026a，黑盒可对标口径）

| 指标 | 数值 |
|---|---|
| **Coverage（找峰率）** | **0.8993** |
| **Precision（精度）** | **0.9921** |
| **F1** | **0.9345** |
| **SR（完全找齐率）** | 0.5383 |
| **官方得分 ½(RPR+F1)** | **0.9169** |

> 官方得分为 0.9169，高于 S-CARD-CMSA（RS-CMSA-ESII + 密度过滤报告）论文公布的 0.818。详见 `REPORT.md` §10。

## 说明

- 全程**黑盒**：不使用真峰位置做优化，生态位半径与密度过滤半径均从候选解自估，与 CEC 官方规则对齐。
- 评估预算固定为 **20000 × dim**（CEC 2026 多模态官方口径）。
- 基准套件 `cec2026/` 移植自 CEC 2026 官方 benchmark（原作者版权归其各自所有）。