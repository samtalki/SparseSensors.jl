### 
# Const-constrained QR optimizer for sensor placement.
###

mutable struct CostQRPivot
    Ψ::AbstractArray
    pivots::AbstractArray
    sensor_costs::AbstractArray #Costs/weights associated with each sensor. Postsive is bad, negative good.
end

function CostQRPivot(qrpivot::QRPivot,sensor_costs::Vector{Float64})
    return CostQRPivot(qrpivot.Ψ,qrpivot.pivots,sensor_costs)
end

function CostQRPivot(Ψ::AbstractArray,pivots::Vector,sensor_costs::Vector{Float64})
    return CostQRPivot(Ψ,pivots,sensor_costs)
end

function fit(cost_qr_pivot::CostQRPivot)
    n,m = size(cost_qr_pivot.Ψ)

    if cost_qr_pivot.sensor_costs !== n
        throw(DomainError(cost_qr_pivot.sensor_costs,"The specified sensor costs are inconsistent with the basis dimensions"))
    end
    #Initialize helper variables
    R,p,k = make_helper_variables(cost_qr_pivot.Ψ)

    for j in 1:k
        u,i_piv = qr_reflector(R[j:n,j:m],cost_qr_pivot.sensor_costs)
        #Track column pivots
        i_piv += j 
        p[[j,i_piv]] = p[[i_piv,j]]
        #Switch column
        R[:,[j,i_piv]] = R[:,[i_piv,j]] 
        #Apply reflector
        R[j:n,j:m] -= outer(u,dot(u,R[j:n,j:m]))
        R[j+1:n,j] = 0
    end
    cost_qr_pivot.pivots = p
end

function make_helper_variables(Ψ)
    """Constructs the reflector helper variables"""
    R = conj(Ψ)
    p = [i for i in 1:n]
    k = min(m,n)
    return R,p,k
end

function qr_reflector(r,costs)
    dlens = sqrt(sum(abs(eachrow(r))^2)) #Norm of each column
    i_piv = argmax(dlens-costs) # choose pivot
    dlen = dlens[i_piv]

    if dlen > 0 
        u = r[:,i_piv]/dlen
        u[1] += sign(u[1])+(u[1]==0)
        u/=sqrt(abs(u[1]))
    else
        u = r[:,i_piv]
        u[1] = sqrt(2)
    end
    return u,i_piv
end
