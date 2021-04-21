using LinearAlgebra, GenericSVD
@testset "Opt common" begin

function gen_test_setting()
    coeffs = [4.0-0.5im, 0.1im, 2.0-0.2im]
    discr = collect(0:9.0)
    n = length(discr)
    A = hcat(ones(n), discr, discr.^2)
    b = coeffs[1] .+ coeffs[2]*discr + coeffs[3]*discr.^2
    return (A,b,coeffs,discr)
end


    # Testing generic complex version
    for linlsqr = [:backslash, :svd, :nrmeq]
        (A,b,coeffs,discr) = gen_test_setting()
        x = solve_linlsqr!(A, b, linlsqr, 0.0)
        @test x ≈ coeffs
        @test !isreal(x)
    end
    (A,b,coeffs,discr) = gen_test_setting()
    x = solve_linlsqr!(A, b, :svd, 1e-1) #testing droptolerance, high = garbage result
    @test norm(x-[4.0, 0.0, 2.0]) > 1


    # Testing to force real solution
    for linlsqr = [:real_backslash, :real_svd, :real_nrmeq]
        (A,b,coeffs,discr) = gen_test_setting()
        A = complex(A)
        x = solve_linlsqr!(A, b, linlsqr, 0.0)
        @test x ≈ real(coeffs)
        @test isreal(x)
    end
    (A,b,coeffs,discr) = gen_test_setting()
    A = complex(A)
    x = solve_linlsqr!(A, b, :svd, 1e-1) #testing droptolerance, high = garbage result
    @test norm(x-[4.0, 0.0, 2.0]) > 1


    # Testing BigFloat arithmetic
    (A,b,coeffs,discr) = gen_test_setting()
    A = big.(A)
    b = big.(b)
    x = solve_linlsqr!(A, b, :svd, 0.0)
    @test x ≈ coeffs

    (A,b,coeffs,discr) = gen_test_setting()
    x = solve_linlsqr!(A, b, :real_svd, 0.0)
    @test x ≈ real(coeffs)
    @test isreal(x)

end
