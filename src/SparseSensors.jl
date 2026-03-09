module SparseSensors
using LinearAlgebra

"""
    AbstractSampler

Abstract supertype for all sensor placement samplers.

Concrete subtypes ([`QRPivot`](@ref), [`CostQRPivot`](@ref)) must provide a
`pivots::Vector{Int}` field and a [`fit`](@ref) method.
"""
abstract type AbstractSampler end

include("qr.jl")
include("ccqr.jl")
include("basis.jl")
include("reconstruct.jl")

export AbstractSampler
export QRPivot, CostQRPivot, fit, VandermondeBasis, SVDBasis, Basis
export get_sensors, reconstruct

"""
    get_sensors(sampler::AbstractSampler) -> Vector{Int}

Return all ranked sensor locations from a fitted [`AbstractSampler`](@ref).
"""
function get_sensors(sampler::AbstractSampler)
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
    get_sensors(sampler::AbstractSampler, n_sensors) -> Vector{Int}

Return the top `n_sensors` locations from a fitted [`AbstractSampler`](@ref).
"""
function get_sensors(sampler::AbstractSampler, n_sensors)
    return get_sensors(sampler.pivots, n_sensors)
end
end
