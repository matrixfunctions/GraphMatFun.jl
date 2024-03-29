@testset "Opt common" begin
    function graph_test_setting(d = 0)
        coeffs = [4.0 - 0.5im, 0.1im, 2.0 - 0.2im]
        discr = collect(0:9.0) .+ d
        n = length(discr)
        A = hcat(ones(n), discr, discr .^ 2)
        b = coeffs[1] .+ coeffs[2] * discr + coeffs[3] * discr .^ 2
        return (A, b, coeffs, discr)
    end

    # Testing generic complex version
    for linlsqr in [:backslash, :svd, :nrmeq]
        (A, b, coeffs, discr) = graph_test_setting()
        x = solve_linlsqr!(A, b, linlsqr, 0.0)
        @test x ≈ coeffs
        @test !isreal(x)
    end
    (A, b, coeffs, discr) = graph_test_setting()
    #testing droptolerance, high = garbage result
    x = solve_linlsqr!(A, b, :svd, 1e-1)
    @test norm(x - [4.0, 0.0, 2.0]) > 1

    # Testing to force real solution
    for linlsqr in [:real_backslash, :real_svd, :real_nrmeq]
        (A, b, coeffs, discr) = graph_test_setting()
        A = complex(A)
        x = solve_linlsqr!(A, b, linlsqr, 0.0)
        @test x ≈ real(coeffs)
        @test isreal(x)

        (A, b, coeffs, discr) = graph_test_setting(0.1im)
        A = complex(A)
        x = solve_linlsqr!(A, b, linlsqr, 0.0)
        @test isreal(x)
        xx = real(A \ b)
        @test norm(A * x - b) < norm(A * xx - b)
    end
    (A, b, coeffs, discr) = graph_test_setting()
    A = complex(A)
    #testing droptolerance, high = garbage result
    x = solve_linlsqr!(A, b, :real_svd, 1e-1)
    @test norm(x - [4.0, 0.0, 2.0]) > 1

    # Testing BigFloat arithmetic
    (A, b, coeffs, discr) = graph_test_setting()
    A = big.(A)
    b = big.(b)
    x = solve_linlsqr!(A, b, :svd, 0.0)
    @test x ≈ coeffs

    (A, b, coeffs, discr) = graph_test_setting()
    A = big.(A)
    b = big.(b)
    x = solve_linlsqr!(A, b, :real_svd, 0.0)
    @test x ≈ real(coeffs)
    @test isreal(x)
end
