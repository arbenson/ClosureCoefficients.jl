using ClosureCoefficients
using SparseArrays
using StatsBase
using Test

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
    @test ret["oo_o"].global_clcf ≈ 1.0
    @test ret["oo_o"].avg_clcf    ≈ 1.0
    @test ret["oo_o"].local_clcfs ≈ [1.0, 0.0, 0.0]
    @test ret["oo_o"].wedges    == [1, 0, 0]
    @test ret["oo_o"].triangles == [1, 0, 0]

    @test ret["oo_i"].global_clcf ≈ 0.0
    @test ret["oo_i"].avg_clcf    ≈ 0.0
    @test ret["oo_i"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["oo_i"].wedges    == [1, 0, 0]
    @test ret["oo_i"].triangles == [0, 0, 0]

    @test ret["oi_o"].global_clcf ≈ 0.5
    @test ret["oi_o"].avg_clcf    ≈ 0.5
    @test ret["oi_o"].local_clcfs ≈ [1.0, 0.0, 0.0]
    @test ret["oi_o"].wedges    == [1, 1, 0]
    @test ret["oi_o"].triangles == [1, 0, 0]

    @test ret["oi_i"].global_clcf ≈ 0.5
    @test ret["oi_i"].avg_clcf    ≈ 0.5
    @test ret["oi_i"].local_clcfs ≈ [0.0, 1.0, 0.0]
    @test ret["oi_i"].wedges    == [1, 1, 0]
    @test ret["oi_i"].triangles == [0, 1, 0]

    @test ret["io_o"].global_clcf ≈ 0.5
    @test ret["io_o"].avg_clcf    ≈ 0.5
    @test ret["io_o"].local_clcfs ≈ [0.0, 1.0, 0.0]
    @test ret["io_o"].wedges    == [0, 1, 1]
    @test ret["io_o"].triangles == [0, 1, 0]

    @test ret["io_i"].global_clcf ≈ 0.5
    @test ret["io_i"].avg_clcf    ≈ 0.5
    @test ret["io_i"].local_clcfs ≈ [0.0, 0.0, 1.0]
    @test ret["io_i"].wedges    == [0, 1, 1]
    @test ret["io_i"].triangles == [0, 0, 1]

    @test ret["ii_o"].global_clcf ≈ 0.0
    @test ret["ii_o"].avg_clcf    ≈ 0.0
    @test ret["ii_o"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["ii_o"].wedges    == [0, 0, 1]
    @test ret["ii_o"].triangles == [0, 0, 0]

    @test ret["ii_i"].global_clcf ≈ 1.0
    @test ret["ii_i"].avg_clcf    ≈ 1.0
    @test ret["ii_i"].local_clcfs ≈ [0.0, 0.0, 1.0]
    @test ret["ii_i"].wedges    == [0, 0, 1]
    @test ret["ii_i"].triangles == [0, 0, 1]
end

function dir_test3()
    # Cycle
    I = [1, 2, 3]
    J = [2, 3, 1]
    A = sparse(I, J, 1, 3, 3)

    ret = dir_clcfs(A)
    @test ret["oo_o"].global_clcf ≈ 0.0
    @test ret["oo_o"].avg_clcf    ≈ 0.0
    @test ret["oo_o"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["oo_o"].wedges    == [1, 1, 1]
    @test ret["oo_o"].triangles == [0, 0, 0]

    @test ret["oo_i"].global_clcf ≈ 1.0
    @test ret["oo_i"].avg_clcf    ≈ 1.0
    @test ret["oo_i"].local_clcfs ≈ [1.0, 1.0, 1.0]
    @test ret["oo_i"].wedges    == [1, 1, 1]
    @test ret["oo_i"].triangles == [1, 1, 1]

    @test ret["oi_o"].global_clcf ≈ 0.0
    @test ret["oi_o"].avg_clcf    ≈ 0.0
    @test ret["oi_o"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["oi_o"].wedges    == [0, 0, 0]
    @test ret["oi_o"].triangles == [0, 0, 0]

    @test ret["oi_i"].global_clcf ≈ 0.0
    @test ret["oi_i"].avg_clcf    ≈ 0.0
    @test ret["oi_i"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["oi_i"].wedges    == [0, 0, 0]
    @test ret["oi_i"].triangles == [0, 0, 0]

    @test ret["io_o"].global_clcf ≈ 0.0
    @test ret["io_o"].avg_clcf    ≈ 0.0
    @test ret["io_o"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["io_o"].wedges    == [0, 0, 0]
    @test ret["io_o"].triangles == [0, 0, 0]

    @test ret["io_i"].global_clcf ≈ 0.0
    @test ret["io_i"].avg_clcf    ≈ 0.0
    @test ret["io_i"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["io_i"].wedges    == [0, 0, 0]
    @test ret["io_i"].triangles == [0, 0, 0]

    @test ret["ii_o"].global_clcf ≈ 1.0
    @test ret["ii_o"].avg_clcf    ≈ 1.0
    @test ret["ii_o"].local_clcfs ≈ [1.0, 1.0, 1.0]
    @test ret["ii_o"].wedges    == [1, 1, 1]
    @test ret["ii_o"].triangles == [1, 1, 1]

    @test ret["ii_i"].global_clcf ≈ 0.0
    @test ret["ii_i"].avg_clcf    ≈ 0.0
    @test ret["ii_i"].local_clcfs ≈ [0.0, 0.0, 0.0]
    @test ret["ii_i"].wedges    == [1, 1, 1]
    @test ret["ii_i"].triangles == [0, 0, 0]
end

function data_test()
    A1 = load_example_data("FW-Florida.txt", symm=false)
    dc = dir_clcfs(A1)
    # TODO(arb): check these numbers more precisely with data from Hao
    @test ≈(dc["ii_i"].avg_clcf, 0.21, atol=0.03)
    @test ≈(dc["oo_i"].avg_clcf, 0.05, atol=0.03)
    A2 = load_example_data("arxiv-AstroPh.txt", symm=true)
    uc = undir_clcfs(A2)
    @test ≈(uc.avg_clcf, 0.250, atol=0.001)
    @test ≈(uc.global_clcf, 0.318, atol=0.001)
end

function run_all()
    undir_test1()
    dir_test1()
    dir_test2()
    dir_test3()
    data_test()
end
run_all()
