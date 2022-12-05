# SEM4CG

Computational algorithm for predicting electric currents from implicit-solvent trajectories of C-alpha/3SPN-based coarse-grained molecular simulations.
SEM4CG computes electric currents by two steps.

1. Discretize 3-D space into cubic grids and calculate a conductivity for each point. Conductivity is based on the function of the distance between a point and closest bead of molecules in the system.
2. Solve partial differential equations (PDEs).

## Usage

```
python3 /path-to-sem/src/main.py --toml meta-file --out output-file -v
```

## Parameters for conductivity

Conductivity functions for amino-acid residues are computed from all-atom MD simulations of a C-alpha-fixed amino acid in constant electrostatic field (100 mV/nm).
For nucleotides, entire heavy atoms of 10 base-pair dsDNA with random sequence are fixed in constant electrostatic field (100 mV/nm).
