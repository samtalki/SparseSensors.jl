# SparseSensors.jl

Data-driven sparse sensor placement for reconstruction, based on
[Manohar et al. (2018)](https://doi.org/10.1109/MCS.2018.2810460).
This package provides a Julia implementation inspired by
[pysensors](https://github.com/dynamicslab/pysensors)
([de Silva et al., 2021](https://doi.org/10.21105/joss.02828)).

## Overview

SparseSensors.jl finds optimal sensor locations from a basis (SVD/POD modes or
Vandermonde polynomials) using column-pivoted QR factorization, then reconstructs
full-state fields from the sparse measurements.

## Quick Start

```julia
using SparseSensors, LinearAlgebra

# Build an SVD basis from snapshot data
X = load_your_data()          # (n_spatial, n_snapshots)
basis = SVDBasis(X, 10)       # keep 10 modes

# Find optimal sensor locations
qr_pivot = QRPivot(basis.Ψ')
fit(qr_pivot)
sensors = get_sensors(qr_pivot, 10)

# Reconstruct a new observation from sparse measurements
y = new_observation[sensors]
x̂ = reconstruct(basis, sensors, y)
```

## API Reference

### Basis Construction

```@docs
Basis
SVDBasis
VandermondeBasis
```

### Sensor Placement

```@docs
AbstractSampler
QRPivot
CostQRPivot
fit
get_sensors
```

### Reconstruction

```@docs
reconstruct
```

### Internal

```@docs
SparseSensors.make_helper_variables
SparseSensors.qr_reflector
```
