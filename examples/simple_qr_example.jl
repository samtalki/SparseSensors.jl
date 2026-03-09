using SparseSensors
using LinearAlgebra

# Setup the experiment
r = 11  # Number of basis modes
n = 1000
x = collect(0.0:1/n:1.0)
basis = VandermondeBasis(x, r)

# Setup QR pivot sensor placement algorithm
qr_pivot = QRPivot(basis.Ψ')
fit(qr_pivot)

# Select the top 15 sensor locations
f = abs.(x .^ 2 .- 0.5)
selected_sensors = get_sensors(qr_pivot, 15)
x_sensed = x[selected_sensors]
y_sensed = f[selected_sensors]

println("Top 15 sensor locations (indices): ", selected_sensors)
println("Sensor x-coordinates: ", round.(x_sensed, digits=4))
println("Sensor measurements:  ", round.(y_sensed, digits=4))

# Reconstruct the signal from sensor measurements
x̂ = reconstruct(basis, selected_sensors, f[selected_sensors])
max_error = maximum(abs.(x̂ .- f))
println("Max reconstruction error: ", max_error)
