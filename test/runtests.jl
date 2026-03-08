using SparseSensors
using LinearAlgebra
using Test

@testset "SparseSensors.jl" begin

    @testset "VandermondeBasis" begin
        x = collect(0.0:0.1:1.0)
        r = 3
        basis = VandermondeBasis(x, r)
        @test size(basis.Ψ) == (length(x), r)
        @test all(basis.Ψ[:, 1] .≈ 1.0)
        @test basis.Ψ[:, 2] ≈ x
        @test basis.Ψ[:, 3] ≈ x .^ 2
    end

    @testset "SVDBasis" begin
        n, m, r = 50, 30, 5
        # Build a rank-r matrix so the first r singular vectors are well-defined
        A = randn(n, r)
        B = randn(r, m)
        X = A * B

        basis = SVDBasis(X, r)
        @test size(basis.Ψ) == (n, r)
        @test basis.r == r

        # Columns should be orthonormal
        @test basis.Ψ' * basis.Ψ ≈ I(r) atol=1e-10

        # The r modes should capture (nearly) all variance of the rank-r matrix
        U_full, S, _ = svd(X)
        captured = sum(S[1:r].^2) / sum(S.^2)
        @test captured ≈ 1.0 atol=1e-10
    end

    @testset "QRPivot" begin
        x = collect(0.0:0.01:1.0)
        r = 5
        basis = VandermondeBasis(x, r)
        Ψ = copy(transpose(basis.Ψ))
        qr_pivot = QRPivot(Ψ)
        fit(qr_pivot)
        pivots = qr_pivot.pivots
        n = size(Ψ, 2)
        @test length(pivots) == n
        @test sort(pivots) == collect(1:n)
    end

    @testset "get_sensors" begin
        x = collect(0.0:0.01:1.0)
        r = 5
        basis = VandermondeBasis(x, r)
        Ψ = copy(transpose(basis.Ψ))
        qr_pivot = QRPivot(Ψ)
        fit(qr_pivot)

        all_sensors = get_sensors(qr_pivot)
        @test length(all_sensors) == size(Ψ, 2)

        top5 = get_sensors(qr_pivot, 5)
        @test length(top5) == 5

        top3_from_pivots = get_sensors(qr_pivot.pivots, 3)
        @test length(top3_from_pivots) == 3
        @test top3_from_pivots == qr_pivot.pivots[1:3]
    end

    @testset "CostQRPivot" begin
        x = collect(0.0:0.01:1.0)
        r = 5
        basis = VandermondeBasis(x, r)
        Ψ = copy(transpose(basis.Ψ))
        n_sensors = size(Ψ, 2)

        # Uniform zero costs → same result as unconstrained QRPivot
        costs = zeros(n_sensors)
        qr_pivot = QRPivot(Ψ)
        fit(qr_pivot)
        ccqr = CostQRPivot(qr_pivot, costs)
        fit(ccqr)
        @test length(ccqr.pivots) == n_sensors
        @test sort(ccqr.pivots) == collect(1:n_sensors)
        @test ccqr.pivots == qr_pivot.pivots

        # High cost on sensor 1 pushes it out of first pivot
        costs_high = zeros(n_sensors)
        costs_high[1] = 1e10
        ccqr2 = CostQRPivot(Ψ, Int[], costs_high)
        fit(ccqr2)
        @test ccqr2.pivots[1] != 1

        # Making one sensor cheap (all others expensive) forces it to be selected first
        costs_cheap = fill(1e10, n_sensors)
        cheap_idx = 50
        costs_cheap[cheap_idx] = 0.0
        ccqr3 = CostQRPivot(Ψ, Int[], costs_cheap)
        fit(ccqr3)
        @test ccqr3.pivots[1] == cheap_idx

        # Wrong-length costs throws DomainError
        @test_throws DomainError fit(CostQRPivot(Ψ, Int[], zeros(3)))
    end

    @testset "reconstruct" begin
        # Reconstruct a known polynomial from sensor measurements
        x = collect(0.0:0.01:1.0)
        n = length(x)
        r = 4

        # True signal: 2 + 3x - x² + 0.5x³
        f_true = 2.0 .+ 3.0 .* x .- x .^ 2 .+ 0.5 .* x .^ 3

        basis = VandermondeBasis(x, r)
        Ψ_T = copy(transpose(basis.Ψ))
        qr_pivot = QRPivot(Ψ_T)
        fit(qr_pivot)

        sensors = get_sensors(qr_pivot, r)
        measurements = f_true[sensors]

        # Reconstruct from explicit sensor indices
        x̂ = reconstruct(basis, sensors, measurements)
        @test length(x̂) == n
        @test x̂ ≈ f_true atol=1e-8

        # Reconstruct via sampler convenience method
        x̂2 = reconstruct(qr_pivot, basis, measurements)
        @test x̂2 ≈ f_true atol=1e-8

        # Unfitted sampler throws ArgumentError
        unfitted = QRPivot(Ψ_T)
        @test_throws ArgumentError reconstruct(unfitted, basis, measurements)
    end

    @testset "SVDBasis + reconstruct" begin
        # Build synthetic data from known modes, place sensors, reconstruct
        n = 100
        t = range(0, 2π, length=n)

        # 3 spatial modes
        mode1 = sin.(t)
        mode2 = cos.(t)
        mode3 = sin.(2t)

        # 50 snapshots with known coefficients
        m = 50
        coeffs = randn(3, m)
        X = hcat(mode1, mode2, mode3) * coeffs  # (n, m)

        basis = SVDBasis(X, 3)
        Ψ_T = copy(transpose(basis.Ψ))
        qr_pivot = QRPivot(Ψ_T)
        fit(qr_pivot)

        # Pick a test snapshot and reconstruct
        test_signal = X[:, 1]
        sensors = get_sensors(qr_pivot, 3)
        measurements = test_signal[sensors]
        x̂ = reconstruct(basis, sensors, measurements)

        @test length(x̂) == n
        @test x̂ ≈ test_signal atol=1e-8
    end
end
