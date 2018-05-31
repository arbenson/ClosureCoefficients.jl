using ClosureCoefficients
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

function undir_test1()
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
        @test ret.global_clcf ≈ global_clcf
        @test ret.avg_clcf    ≈ avg_clcf
        @test ret.local_clcfs ≈ local_clcfs      
        @test ret.wedges     == wedges
        @test ret.triangles  == triangles
    end
end

function dir_test1()
    # Simple all bidirectional edges
    I = [1, 1, 2, 2, 3, 3]
    J = [2, 3, 1, 3, 1, 2]
    A = sparse(I, J, 1, 3, 3)

    ret = dir_clcfs(A)
    for (k, v) in ret
        @test v.global_clcf ≈ 1.0
        @test v.avg_clcf    ≈ 1.0
        @test v.local_clcfs ≈ [1.0, 1.0, 1.0]
        @test v.wedges    == [2, 2, 2]
        @test v.triangles == [2, 2, 2]
    end
end

function dir_test2()
    # Feed forward loop
    I = [1, 1, 2]
    J = [2, 3, 3]
    A = sparse(I, J, 1, 3, 3)

    ret = dir_clcfs(A)
    @test ret["OO_O"].global_clcf ≈ 1.0
    @test ret["OO_O"].avg_clcf    ≈ 1.0
    @test ret["OO_O"].local_clcfs ≈ [1.0, 0.0, 0.0]
    @test ret["OO_O"].wedges    == [1, 0, 0]
    @test ret["OO_O"].triangles == [1, 0, 0]

    @test ret["OO_I"].global_clcf ≈ 0.0
    @test ret["OO_I"].avg_clcf    ≈ 0.0
    @test ret["OO_I"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["OO_I"].wedges    == [1, 0, 0]
    @test ret["OO_I"].triangles == [0, 0, 0]

    @test ret["OI_O"].global_clcf ≈ 0.5
    @test ret["OI_O"].avg_clcf    ≈ 0.5
    @test ret["OI_O"].local_clcfs ≈ [1.0, 0.0, 0.0]
    @test ret["OI_O"].wedges    == [1, 1, 0]
    @test ret["OI_O"].triangles == [1, 0, 0]

    @test ret["OI_I"].global_clcf ≈ 0.5
    @test ret["OI_I"].avg_clcf    ≈ 0.5
    @test ret["OI_I"].local_clcfs ≈ [0.0, 1.0, 0.0]
    @test ret["OI_I"].wedges    == [1, 1, 0]
    @test ret["OI_I"].triangles == [0, 1, 0]

    @test ret["IO_O"].global_clcf ≈ 0.5
    @test ret["IO_O"].avg_clcf    ≈ 0.5
    @test ret["IO_O"].local_clcfs ≈ [0.0, 1.0, 0.0]
    @test ret["IO_O"].wedges    == [0, 1, 1]
    @test ret["IO_O"].triangles == [0, 1, 0]

    @test ret["IO_I"].global_clcf ≈ 0.5
    @test ret["IO_I"].avg_clcf    ≈ 0.5
    @test ret["IO_I"].local_clcfs ≈ [0.0, 0.0, 1.0]
    @test ret["IO_I"].wedges    == [0, 1, 1]
    @test ret["IO_I"].triangles == [0, 0, 1]

    @test ret["II_O"].global_clcf ≈ 0.0
    @test ret["II_O"].avg_clcf    ≈ 0.0
    @test ret["II_O"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["II_O"].wedges    == [0, 0, 1]
    @test ret["II_O"].triangles == [0, 0, 0]

    @test ret["II_I"].global_clcf ≈ 1.0
    @test ret["II_I"].avg_clcf    ≈ 1.0
    @test ret["II_I"].local_clcfs ≈ [0.0, 0.0, 1.0]
    @test ret["II_I"].wedges    == [0, 0, 1]
    @test ret["II_I"].triangles == [0, 0, 1]
end

function dir_test3()
    # Cycle
    I = [1, 2, 3]
    J = [2, 3, 1]
    A = sparse(I, J, 1, 3, 3)

    ret = dir_clcfs(A)
    @test ret["OO_O"].global_clcf ≈ 0.0
    @test ret["OO_O"].avg_clcf    ≈ 0.0
    @test ret["OO_O"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["OO_O"].wedges    == [1, 1, 1]
    @test ret["OO_O"].triangles == [0, 0, 0]

    @test ret["OO_I"].global_clcf ≈ 1.0
    @test ret["OO_I"].avg_clcf    ≈ 1.0
    @test ret["OO_I"].local_clcfs ≈ [1.0, 1.0, 1.0]
    @test ret["OO_I"].wedges    == [1, 1, 1]
    @test ret["OO_I"].triangles == [1, 1, 1]

    @test ret["OI_O"].global_clcf ≈ 0.0
    @test ret["OI_O"].avg_clcf    ≈ 0.0
    @test ret["OI_O"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["OI_O"].wedges    == [0, 0, 0]
    @test ret["OI_O"].triangles == [0, 0, 0]

    @test ret["OI_I"].global_clcf ≈ 0.0
    @test ret["OI_I"].avg_clcf    ≈ 0.0
    @test ret["OI_I"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["OI_I"].wedges    == [0, 0, 0]
    @test ret["OI_I"].triangles == [0, 0, 0]

    @test ret["IO_O"].global_clcf ≈ 0.0
    @test ret["IO_O"].avg_clcf    ≈ 0.0
    @test ret["IO_O"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["IO_O"].wedges    == [0, 0, 0]
    @test ret["IO_O"].triangles == [0, 0, 0]

    @test ret["IO_I"].global_clcf ≈ 0.0
    @test ret["IO_I"].avg_clcf    ≈ 0.0
    @test ret["IO_I"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["IO_I"].wedges    == [0, 0, 0]
    @test ret["IO_I"].triangles == [0, 0, 0]

    @test ret["II_O"].global_clcf ≈ 1.0
    @test ret["II_O"].avg_clcf    ≈ 1.0
    @test ret["II_O"].local_clcfs ≈ [1.0, 1.0, 1.0]
    @test ret["II_O"].wedges    == [1, 1, 1]
    @test ret["II_O"].triangles == [1, 1, 1]

    @test ret["II_I"].global_clcf ≈ 0.0
    @test ret["II_I"].avg_clcf    ≈ 0.0
    @test ret["II_I"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["II_I"].wedges    == [1, 1, 1]
    @test ret["II_I"].triangles == [0, 0, 0]
end

function run_all()
    undir_test1()
    dir_test1()
    dir_test2()
    dir_test3()
end
run_all()
