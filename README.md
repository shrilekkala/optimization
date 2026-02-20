# Optimization (MATLAB): Augmented Lagrangian + Bound-Constrained Subproblems

This repository contains MATLAB implementations from an optimization final project (STAT 310). The core work solves randomly generated nonlinear programs with **equality constraints** and **box constraints** using an **augmented Lagrangian / bound-constrained Lagrangian method**, and evaluates runtime scaling in convex and nonconvex settings.

Final submission: **[Code](./FINAL/)** · **[Report](FINAL/Optimization%20Final%20-%20Augmented%20Lagrangian%20%2B%20Bound-Constrained%20Subproblems.pdf)**

---

## Problem (High Level)
We solve randomly generated nonlinear programs with:
- **Equality constraints**
- **Bound constraints:** -1 <= x <= 1

A parameter \(\gamma\) controls nonconvexity (\(\gamma = 0\) gives the convex quadratic case).

---

## Approach
- Outer loop: **Bound-Constrained Lagrangian / Augmented Lagrangian** iterations
- Inner loop: solve bound-constrained subproblems using **projected methods** (gradient projection + projected CG)
- Includes utilities for:
  - Cauchy point computation
  - projection onto bounds
  - SOSC checking (for the nonconvex setting)

---

## Repository Structure

### `FINAL/` (main project code + experiments)
**Entry point**
- `final.m` — required solver entry:
  - `[x,obj] = final(m,n,s,p)`

**Core solver components**
- `BCLMethod.m`
- `augmentedLagrangian.m`
- `solveSubproblem.m`
- `gradLagrangian.m`
- `computeHessian.m`

**Bound-constrained / projected optimization**
- `gradient_projection_method.m`
- `computeCauchyPoint.m`
- `projected_CG.m`
- `projectOntoBounds.m`

**Analysis / experiments**
- `final_p2.m` — convex scaling experiment (γ = 0)
- `final_p3.m` — nonconvex scaling experiment (γ > 0)
- `final_p4.m` — SOSC frequency experiment
- `checkSOSC.m`

**Outputs**
- `figures/` — generated plots

---

## How to Run

### 1) Open MATLAB and add `FINAL/` to path
From the repo root:
```
addpath("FINAL")
```

### 2) Run the solver once (example)
```
[x,obj] = final(5,10,6,100);
```

### 3) Reproduce the experiments/figures
```
final_p2   % convex scaling (gamma = 0)
final_p3   % nonconvex scaling (gamma > 0)
final_p4   % SOSC frequency check
```

---

## Notes
- The scripts are written for reproducibility of the final project experiments and plots.
- If the solver hits the iteration cap p before meeting tolerances, it returns the last iterate.
