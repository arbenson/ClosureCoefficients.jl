using ClosureCoefficients
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

function basic_undir_test1()
    I = [1, 1, 2, 3, 3]
    J = [2, 3, 3, 4, 5]
    A = sparse(I, J, 1, 5, 5)

    global_clcf = (2 + 2 + 2 + 0 + 0) / (4 + 4 + 2 + 3 + 3)
    avg_clcf    = mean([1.0 / 2, 1.0 / 2, 1.0, 0.0, 0.0])
    local_clcfs =      [1.0 / 2, 1.0 / 2, 1.0, 0.0, 0.0]
    wedges      = [4, 4, 2, 3, 3]
    triangles   = [1, 1, 1, 0, 0]    

    C = copy(A)
    C[1, 1] = 1
    for B in [max.(A, A'), A, C, max.(C, C')]
        ret = undir_clcfs(B)
        @test ret.global_clcf  ≈ global_clcf
        @test ret.avg_clcf     ≈ avg_clcf
        @test ret.local_clcfs  ≈ local_clcfs      
        @test ret.wedges == wedges
        @test ret.triangles  == triangles
    end
end

function run()
    basic_undir_test1()
end
run()
