###
# Standard QR Pivot Optimizer
###

mutable struct QRPivot
	basis_matrix::Matrix
	pivots::AbstractArray #Ranked list of sensor locations
	fit::Function
end

function QRPivot(Ψ::Matrix,pivots)
	return QRPivot(Ψ,pivots)
end

function fit(Ψ::Matrix;optimizer_kwargs...)
	"""
	Fits the QRPivot sensor placement.
	Ψ: Matrix
		Basis matrix
	optimizer_kwargs: dictionary
		Optional settings for optimizer
	"""
	pivots = sensor_placement(Ψ)
	res = QRPivot(Ψ,pivots)
	return res
end

function get_sensors(sparse_sampler::QRPivot)
	return sprase_sampler.pivots
end

function sensor_placement(Ψ)
	# --> Compute the QRPivot w/ column pivoting decomposition of Ψ.
	_, _, p = qr(transpose(Ψ), Val(true))
	return p[1:size(Ψ, 2)]
end


