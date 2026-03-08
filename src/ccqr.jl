"""
    CostQRPivot(Ψ, pivots, sensor_costs)

Cost-constrained QR sensor placement optimizer.

Extends QR-pivoted sensor placement by penalizing expensive sensor locations.
Sensors with high cost are pushed later in the pivot ordering.

# Fields
- `Ψ::AbstractArray`: Basis matrix (modes × sensors).
- `pivots::Vector{Int}`: Ranked sensor locations (populated by [`fit`](@ref)).
- `sensor_costs::AbstractVector`: Per-sensor cost vector (length must equal number of columns).

See also: [`QRPivot`](@ref), [`fit`](@ref)
"""
mutable struct CostQRPivot <: AbstractSampler
    Ψ::AbstractArray
    pivots::Vector{Int}
    sensor_costs::AbstractVector
end

function CostQRPivot(qrpivot::QRPivot, sensor_costs::AbstractVector)
    return CostQRPivot(qrpivot.Ψ, Int[], sensor_costs)
end

"""
    fit(cost_qr_pivot::CostQRPivot) -> CostQRPivot

Compute cost-constrained sensor locations. Works like [`fit(::QRPivot)`](@ref) but
subtracts `sensor_costs` from column norms at each pivot step, preferring cheaper sensors
when information content is similar.

Throws `DomainError` if the length of `sensor_costs` does not match the number of columns.
"""
function fit(cost_qr_pivot::CostQRPivot)
    n, m = size(cost_qr_pivot.Ψ)

    if length(cost_qr_pivot.sensor_costs) != m
        throw(DomainError(cost_qr_pivot.sensor_costs, "Sensor costs length must match number of columns (sensors)"))
    end

    R, p, k = make_helper_variables(cost_qr_pivot.Ψ)

    for j in 1:k
        Rsub = @view R[j:n, j:m]
        u, i_piv = qr_reflector(Rsub, cost_qr_pivot.sensor_costs[p[j:m]])
        i_piv += j - 1
        p[[j, i_piv]] = p[[i_piv, j]]
        R[:, [j, i_piv]] = R[:, [i_piv, j]]
        Rsub .-= u * (u' * Rsub)
        R[j+1:n, j] .= 0
    end
    cost_qr_pivot.pivots = p
    return cost_qr_pivot
end

"""
    make_helper_variables(Ψ) -> (R, p, k)

Initialize working variables for cost-constrained QR factorization.

Returns the conjugate copy `R`, the identity permutation `p`, and the pivot count `k`.
"""
function make_helper_variables(Ψ)
    R = float.(conj(Ψ))
    n, m = size(Ψ)
    p = collect(1:m)
    k = min(m, n)
    return R, p, k
end

"""
    qr_reflector(r, costs) -> (u, i_piv)

Compute a Householder reflector for cost-constrained column pivoting.

Selects the pivot column that maximises `‖column‖ − cost`, then returns the
Householder vector `u` and the pivot index `i_piv`.
"""
function qr_reflector(r, costs)
    dlens = [norm(c) for c in eachcol(r)]
    i_piv = argmax(dlens - costs)
    dlen = dlens[i_piv]

    if dlen > 0
        u = r[:, i_piv] / dlen
        u[1] += sign(u[1]) + (u[1] == 0)
        u /= sqrt(abs(u[1]))
    else
        u = r[:, i_piv]
        u[1] = sqrt(2)
    end
    return u, i_piv
end
