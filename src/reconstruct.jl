"""
    reconstruct(basis::Basis, sensor_indices::AbstractVector{Int}, measurements::AbstractVector)

Reconstruct the full state from sparse sensor measurements.

Solves the linear system `Ψ[sensor_indices, :] * a = measurements` for the
mode coefficients `a`, then returns the full-field estimate `Ψ * a`.

# Arguments
- `basis::Basis`: The basis used for sensor placement.
- `sensor_indices::AbstractVector{Int}`: Indices of the sensor locations.
- `measurements::AbstractVector`: Measured values at those locations.

# Returns
- `Vector`: Reconstructed full-state vector of length `size(basis.Ψ, 1)`.

# Examples
```julia
basis = VandermondeBasis(0.0:0.01:1.0, 5)
sensors = [1, 26, 51, 76, 101]
y = f.(x[sensors])          # measure at sensor locations
x̂ = reconstruct(basis, sensors, y)
```
"""
function reconstruct(basis::Basis, sensor_indices::AbstractVector{Int}, measurements::AbstractVector)
    Θ = basis.Ψ[sensor_indices, :]
    a = Θ \ measurements
    return basis.Ψ * a
end

"""
    reconstruct(sampler, basis::Basis, measurements::AbstractVector)

Convenience method that extracts sensor indices from `sampler.pivots` (a fitted
[`QRPivot`](@ref) or [`CostQRPivot`](@ref)) and reconstructs the full state.

The number of sensors used equals the number of basis modes (`basis.r`).
"""
function reconstruct(sampler, basis::Basis, measurements::AbstractVector)
    sensor_indices = sampler.pivots[1:basis.r]
    return reconstruct(basis, sensor_indices, measurements)
end
