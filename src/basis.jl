"""
    Basis(Ψ, r)

A truncated basis representation for sparse sensor placement.

# Fields
- `Ψ::Matrix`: Basis matrix of size `(n, r)` where `n` is the number of spatial locations
  and `r` is the number of retained modes.
- `r::Int`: Number of retained modes (columns of `Ψ`).
"""
struct Basis
    Ψ::Matrix
    r::Int
end

"""
    VandermondeBasis(x, r::Int) -> Basis

Construct a [`Basis`](@ref) from a Vandermonde matrix. Column `j` is `x.^(j-1)`.

# Arguments
- `x`: Vector of evaluation points.
- `r::Int`: Number of basis modes (polynomial degree + 1).

# Examples
```julia
basis = VandermondeBasis(0.0:0.1:1.0, 3)
basis.Ψ[:, 1]  # all ones
basis.Ψ[:, 2]  # x values
basis.Ψ[:, 3]  # x² values
```
"""
function VandermondeBasis(x, r::Int)
    Ψ = zeros((length(x), r))
    for j in 1:r
        Ψ[:, j] = make_vander_column(j, x)
    end
    return Basis(Ψ, r)
end

"""
    SVDBasis(X::AbstractMatrix, r::Int) -> Basis

Construct a [`Basis`](@ref) from the first `r` left singular vectors of a data matrix `X`.

This is the primary basis in the sparse sensor placement literature
(Manohar et al., 2018). Each column of `X` is a snapshot; the returned basis
captures the dominant `r` spatial modes.

# Arguments
- `X::AbstractMatrix`: Data matrix of size `(n, m)` where `n` is the number of spatial
  locations and `m` is the number of snapshots.
- `r::Int`: Number of modes to retain (`r ≤ min(n, m)`).

# Examples
```julia
X = randn(100, 50)   # 100 spatial points, 50 snapshots
basis = SVDBasis(X, 5)
size(basis.Ψ)         # (100, 5)
```
"""
function SVDBasis(X::AbstractMatrix, r::Int)
    n, m = size(X)
    r ≤ min(n, m) || throw(ArgumentError("r=$r exceeds min(size(X)) = $(min(n, m))"))
    U, _, _ = svd(X)
    return Basis(U[:, 1:r], r)
end

function make_vander_column(k, col_k)
    ψ = [x^(k-1) for x in col_k]
    return ψ
end
