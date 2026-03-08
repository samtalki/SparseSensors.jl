module SparseSensors
using LinearAlgebra
include("qr.jl")
include("ccqr.jl")
include("basis.jl")
include("reconstruct.jl")

export QRPivot, CostQRPivot, fit, VandermondeBasis, SVDBasis, Basis
export get_sensors, reconstruct

"""
    get_sensors(sampler) -> Vector{Int}

Return all ranked sensor locations from a fitted [`QRPivot`](@ref) or [`CostQRPivot`](@ref).
"""
function get_sensors(sampler)
    return sampler.pivots
end

"""
    get_sensors(pivots::AbstractArray, n_sensors::Int) -> Vector{Int}

Return the top `n_sensors` locations from a pivot vector.
"""
function get_sensors(pivots::AbstractArray, n_sensors::Int)
    return pivots[1:n_sensors]
end

"""
    get_sensors(sampler, n_sensors) -> Vector{Int}

Return the top `n_sensors` locations from a fitted sampler.
"""
function get_sensors(sampler, n_sensors)
    return get_sensors(sampler.pivots, n_sensors)
end
end
