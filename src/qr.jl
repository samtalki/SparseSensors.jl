"""
    QRPivot(Ψ)

Optimal sensor placement via QR factorization with column pivoting.

# Fields
- `Ψ::AbstractArray`: Basis matrix (typically modes × sensors)
- `pivots::Vector{Int}`: Ranked sensor locations (populated by [`fit`](@ref))
"""
mutable struct QRPivot
	Ψ::AbstractArray
	pivots::Vector{Int}
end

function QRPivot(Ψ)
	return QRPivot(Ψ, Int[])
end

"""
    fit(qr_pivot::QRPivot) -> QRPivot

Compute optimal sensor locations for the given [`QRPivot`](@ref) using column-pivoted QR
factorization. The resulting pivot indices are stored in `qr_pivot.pivots`, ordered from
most to least informative.
"""
function fit(qr_pivot::QRPivot)
	qr_pivot.pivots = sensor_placement(qr_pivot.Ψ)
	return qr_pivot
end

function sensor_placement(Ψ)
	_, _, p = qr(conj(Ψ), ColumnNorm())
	return p[1:size(Ψ, 2)]
end
