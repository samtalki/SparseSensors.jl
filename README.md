# SparseSensors.jl

[![CI](https://github.com/samtalki/SparseSensors.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/samtalki/SparseSensors.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://samtalki.github.io/SparseSensors.jl/dev/)

Data-driven sparse sensor placement for reconstruction, implementing the core algorithms from Manohar, *et al.*, "Data-Driven Sparse Sensor Placement for Reconstruction" in Julia. This is a hobbyist attempt to provide a faithful implementation of the fantastic [pysensors](https://github.com/dynamicslab/pysensors) library in Julia.

All collaborations and contributions are welcome.

## Installation
To install, use Pkg. From the Julia REPL, press `]` to enter Pkg-mode and run
```julia
pkg> add SparseSensors
```

## Example

```julia
using SparseSensors
using LinearAlgebra

# Build an SVD basis from snapshot data
X = randn(200, 50)              # 200 spatial points, 50 snapshots
basis = SVDBasis(X, 10)         # keep 10 modes

# Find optimal sensor locations
qr_pivot = QRPivot(basis.Ψ')
fit(qr_pivot)
sensors = get_sensors(qr_pivot, 10)

# Reconstruct a new observation from sparse measurements
x_new = X[:, 1]                 # test signal
y = x_new[sensors]              # measure at sensor locations
x̂ = reconstruct(basis, sensors, y)
```

A Vandermonde (polynomial) basis is also available:

```julia
x = collect(0.0:0.01:1.0)
basis = VandermondeBasis(x, 5)
qr_pivot = QRPivot(basis.Ψ')
fit(qr_pivot)
sensors = get_sensors(qr_pivot, 5)

f = @. 2.0 + 3.0*x - x^2       # true signal
y = f[sensors]
x̂ = reconstruct(basis, sensors, y)
```

## References

- Manohar, Krithika, Bingni W. Brunton, J. Nathan Kutz, and Steven L. Brunton.
  "Data-driven sparse sensor placement for reconstruction: Demonstrating the
  benefits of exploiting known patterns."
  IEEE Control Systems Magazine 38, no. 3 (2018): 63-86.
  [DOI](https://doi.org/10.1109/MCS.2018.2810460)

- Clark, Emily, Travis Askham, Steven L. Brunton, and J. Nathan Kutz.
  "Greedy sensor placement with cost constraints." IEEE Sensors Journal 19, no. 7
  (2018): 2642-2656.
  [DOI](https://doi.org/10.1109/JSEN.2018.2887044)

- de Silva, Brian M., Krithika Manohar, Emily Clark, Bingni W. Brunton,
  Steven L. Brunton, J. Nathan Kutz.
  "PySensors: A Python package for sparse sensor placement."
  arXiv preprint arXiv:2102.13476 (2021). [arXiv](https://arxiv.org/abs/2102.13476)
