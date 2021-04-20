using LinearAlgebra
@testset "Opt common" begin

    coeffs = [4.0, 0.0, 2.0]
    discr = collect(0:9.0)
    n = length(discr)
    A = hcat(ones(n), discr, discr.^2)
    b = coeffs[1] .+ coeffs[2]*discr + coeffs[3]*discr.^2

    for linsqr = [:backslash, :nrmeq, :svd]
        x = solve_linlsqr(A, b, linsqr, 0.0)
        @test x ≈ coeffs
    end
    x = solve_linlsqr(A, b, :svd, 1e-1)
    @test norm(x-[4.0, 0.0, 2.0]) > 1

    coeffs = [4.0-0.5im, 0.1im, 2.0-0.2im]
    A = complex(A)
    b = coeffs[1] .+ coeffs[2]*discr + coeffs[3]*discr.^2
    for linsqr = [:real_backslash, :real_nrmeq]
        x = solve_linlsqr(A, b, linsqr, 0.0)
        @test x ≈ real(coeffs)
        @test isreal(x)
    end

end
