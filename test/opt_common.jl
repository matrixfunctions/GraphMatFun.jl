using LinearAlgebra, GenericSVD
@testset "Opt common" begin

    coeffs = [4.0, 0.0, 2.0]
    discr = collect(0:9.0)
    n = length(discr)
    A = hcat(ones(n), discr, discr.^2)
    b = coeffs[1] .+ coeffs[2]*discr + coeffs[3]*discr.^2

    for linlsqr = [:backslash, :svd, :nrmeq]
        x = solve_linlsqr(A, b, linlsqr, 0.0)
        @test x ≈ coeffs
    end
    x = solve_linlsqr(A, b, :svd, 1e-1) #testing droptolerance, high = garbage result
    @test norm(x-[4.0, 0.0, 2.0]) > 1

    coeffs = [4.0-0.5im, 0.1im, 2.0-0.2im]
    A = complex(A)
    b = coeffs[1] .+ coeffs[2]*discr + coeffs[3]*discr.^2
    for linlsqr = [:real_backslash, :real_svd, :real_nrmeq]
        x = solve_linlsqr(A, b, linlsqr, 0.0)
        @test x ≈ real(coeffs)
        @test isreal(x)
        @test isa(A,Matrix{ComplexF64})
        @test isa(b,Vector{ComplexF64})
    end
    x = solve_linlsqr(A, b, :svd, 1e-1) #testing droptolerance, high = garbage result
    @test norm(x-[4.0, 0.0, 2.0]) > 1


    A = big.(A)
    b = big.(b)
    x = solve_linlsqr(copy(A), b, :svd, 0.0)
    @test x ≈ coeffs
    @test isa(A,Matrix{Complex{BigFloat}})
    @test isa(b,Vector{Complex{BigFloat}})

    x = solve_linlsqr(copy(A), b, :real_svd, 0.0)
    @test x ≈ real(coeffs)
    @test isreal(x)
    @test isa(A,Matrix{Complex{BigFloat}})
    @test isa(b,Vector{Complex{BigFloat}})

end
