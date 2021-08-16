module SparseSensors
using LinearAlgebra: Matrix, conj
greet() = print("Hello World!")
using LinearAlgebra
include("qr.jl")
include("ccqr.jl")
include("basis.jl")
include("reconstruct.jl")

export QRPivot,CostQRPivot,fit,VandermondeBasis

export get_sensors

function get_sensors(sampler)
    return sampler.pivots
end

function get_sensors(pivots::AbstractArray,n_sensors::Int)
    return pivots[1:n_sensors]
end

function get_sensors(sampler,n_sensors)
    return sampler.pivots[1:n_sensors]
end
end

