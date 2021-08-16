# SparseSensors.jl

This repository is an implementation of the core QR factorization and cost-constrained QR factorization algorithms for sparse sensor placement described in Manohar, et al., *Data-Driven Sparse Sensor Placcement for Reconstruction*. This is a hobbyist attempt to provide a faithful implementation of the fantastic pysensors library in Julia. 

All collaborations and contributions are welcome.


## Dependencies

LinearAlgebra. Future versions will likely depend on JuMP, Statistics, and DataFrames.

## Example

```
using SparseSensors
using LinearAlgebra
using Gadfly


```

## Todo:

- Implement high level SSPOR and SSPOC interfaces
- Implement high level basis representation wrapper 
- Set up compatibility with JuliaML/MLJ.jl/ScikitLearn.jl

References
------------

-  Manohar, Krithika, Bingni W. Brunton, J. Nathan Kutz, and Steven L. Brunton.
   "Data-driven sparse sensor placement for reconstruction: Demonstrating the
   benefits of exploiting known patterns."
   IEEE Control Systems Magazine 38, no. 3 (2018): 63-86.
   `[DOI] <https://doi.org/10.1109/MCS.2018.2810460>`

-  Clark, Emily, Travis Askham, Steven L. Brunton, and J. Nathan Kutz.
   "Greedy sensor placement with cost constraints." IEEE Sensors Journal 19, no. 7
   (2018): 2642-2656.
   `[DOI] <https://doi.org/10.1109/JSEN.2018.2887044>`

-  de Silva, Brian M., Krithika Manohar, Emily Clark, Bingni W. Brunton,
   Steven L. Brunton, J. Nathan Kutz.
   "PySensors: A Python package for sparse sensor placement."
   arXiv preprint arXiv:2102.13476 (2021). `[arXiv] <https://arxiv.org/abs/2102.13476>`