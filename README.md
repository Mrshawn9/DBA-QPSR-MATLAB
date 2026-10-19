# DBA-QPSR (MATLAB)

Core MATLAB implementation of the algorithm described in the manuscript:

**DBA-QPSR: A Multi-Strategy Enhanced Detective Behavior Algorithm and Its Application to China's Oil Production Forecasting**  
Xuewen Zhang et al.

## Files

- `DBA_QPSR.m`: Core optimizer integrating PIS, DBA with OBL, Soft-rime, and quadratic interpolation (QI).
- `initialization.m`: Population initialization function.
- `demo_dba_qpsr.m`: Minimal runnable example on a 10-dimensional Sphere function.

## Run

Place the three `.m` files in the same folder. Open MATLAB in that folder and run:

```matlab
demo_dba_qpsr
```

You can also call the optimizer directly:

```matlab
rng(42, 'twister');
fobj = @(x) sum(x.^2);
[bestFitness, bestPosition, convergenceCurve, nfe] = ...
    DBA_QPSR(30, 100, -100, 100, 10, fobj);
```

Inputs: population size, maximum iterations, lower bound, upper bound, dimension, and minimization objective function. The population size must be at least 3. Bounds can be scalars or 1-by-D row vectors.

Outputs: best objective value, best solution vector, convergence history, and number of objective-function evaluations.

## Scope

This repository contains the **core DBA-QPSR algorithm** and a basic usage example. It does **not** include the complete CEC2017/CEC2022 experiments, third-party comparison algorithms, reported experimental result tables, or the oil-production forecasting application. The demonstration is not intended to reproduce the paper's reported benchmark or application results.

## Method references

DBA-QPSR integrates previously published search mechanisms:

- Cheng, J.; De Waele, W. Detective Behavior Algorithm (DBA): A New Metaheuristic for Design and Engineering Optimization. *Knowledge-Based Systems* **2026**, **338**, 115434. https://doi.org/10.1016/j.knosys.2026.115434
- Abdel-Basset, M. et al. Parameter estimation of photovoltaic models using an improved marine predators algorithm. *Energy Conversion and Management* **2021**, **227**, 113491. (PIS)
- Su, H. et al. RIME: A physics-based optimization. *Neurocomputing* **2023**, **532**, 183–214. https://doi.org/10.1016/j.neucom.2023.02.010 (Soft-rime)
- Hu, G. et al. An enhanced manta ray foraging optimization algorithm for shape optimization of complex CCG-Ball curves. *Knowledge-Based Systems* **2022**, **240**, 108071. https://doi.org/10.1016/j.knosys.2021.108071 (QI)

See the accompanying manuscript for methodological details and the full bibliography.
