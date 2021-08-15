###
# Standard QR Pivot Optimizer
###

mutable struct QRPivot
	Ψ::AbstractArray #Basis matrix from SVD, RPCA, etc.
	pivots::AbstractArray #Ranked list of sensor locations
end

function QRPivot(Ψ)
	pivots = zeros(max(size(Ψ)))
	return QRPivot(Ψ,pivots)
end

function QRPivot(Ψ,pivots::AbstractArray)
	return QRPivot(Ψ,pivots)
end

function fit(qrpivot::QRPivot;optimizer_kwargs...)
	"""
	Fits the QRPivot sensor placement.
	qrpivot::QRPivot
		QRPivot object.
	optimizer_kwargs: dictionary
		Optional settings for optimizer
	"""
	qr_pivot.pivots = sensor_placement(qr_pivot.Ψ)
end

function get_sensors(sparse_sampler)
	return sparse_sampler.pivots
end

function sensor_placement(Ψ)
	# --> Compute the QRPivot w/ column pivoting decomposition of Ψ.
	_, _, p = qr(transpose(Ψ), Val(true))
	return p[1:size(Ψ, 2)]
end


