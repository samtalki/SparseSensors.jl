struct Basis
    Ψ::Matrix
    r::Int
end

function VandermondeBasis(x,r::Int;verbose::Bool=false)
    """Makes a vandermonde basis given a range and a number of basis modes"""
    Ψ = zeros((length(x),r))
    for (j,col) in enumerate(eachcol(Ψ))
        ψ = make_vander_column(j,x)
        if verbose
            println("ψ","_"*string(j)," size: ",size(ψ))
        end
        Ψ[:,j] = make_vander_column(j,x) #Each col is x^(k-1)
    end
    if verbose
        print("Final Ψ size: ",size(Ψ))
    end
    return Basis(Ψ,r) #Truncated basis matrix
end

function make_vander_column(k,col_k)
    ψ = [x^(k-1) for x in col_k]
    return ψ #Column of basis matrix
end
