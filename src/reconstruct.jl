using Convex,SCS

mutable struct SSPOR
    basis::Basis
    selected_sensors::Vector
    n_sensors::Int
    n_basis_modes::Int
    fit_::Function
    predict_::Function
    optimizer
end

struct Reconstructor
    fit::Function
    predict::Function
end

function fit(model::SSPOR;quiet=true)
    if model.fit_ === fit_lsq
        model.fit_()
    end
end

function predict(model::SSPOR,x::Vector;quiet=true)
    if model.predict_ === predict_lsq
        model.predict_(model.Ψ,x,model.n_basis_modes)
    end

end

function fit_lsq(Ψ,x,y,selected_sensors)
    y_sensed = y[selected_sensors];
    x_sensed = x[selected_sensors];
    Ψ_sensed = Ψ[:,selected_sensors]
    return Ψ_sensed,x_sensed,y_sensed
end

function predict_lsq(Ψ,x,y,n_basis_modes)
    X = Convex.Variable(size(x))    
    prob = minimize(sumsquares(y-Ψ*X))
    solve!(prob,()-> SCS.Optimizer(verbose=false))
    prob.status
    return prob.optval
end

