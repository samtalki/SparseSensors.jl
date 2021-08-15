

struct VanderMondeBasis
    Ψ::AbstractArray
    r::Int
end

struct SVDBasis
    Ψ::AbstractArray
    r::Int
end

function svd_basis(x,r::Int)
    #To be implemented
end

function vandermonde_basis(x,r::Int)
    """Makes a vandermonde basis given a range and a number of basis modes"""
    Ψ = zeros(Float64,length(x),r)
    for (j,col) in enumerate(eachcol(Ψ))
        Ψ[:,j] = make_vander_column(j,x) #Each col is x^(k-1)
    return VanderMondeBasis(Ψ,r) #Truncated basis matrix
end

function make_vander_column(k,col_k)
    ψ = [x^(k-1) for x in col_k]
    return ψ #Column of basis matrix
end
