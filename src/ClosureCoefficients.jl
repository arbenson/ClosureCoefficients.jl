module ClosureCoefficients

nz_row_inds(A::SpIntMat, ind::Int64) = A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

"""
undir_clcfs_data
----------
This is the data type returned by the function that computes undirected closure coefficients.

The field values are
global_clcf::Float64
    The global closure coefficient.
avg_clcf::Float64
    The average closure coefficient (the mean is taken over nodes that head at least one wedge).
local_clcfs::Vector{Float64}
    A vector of local closure coefficients. If a node is not the head of at least one wedge, then the value is 0.
wedge_counts::Vector{Int64}
    Higher-order wedge counts of nodes: ho_wedge_counts[v] is the number of higher-order wedges with node v at its center.
tri_counts::Vector{Int64}
    Vector of triangle counts of nodes: tri_counts[v] is the number of triangles containing v.
"""
immutable undir_clcfs_data
    global_clcf::Float64
    avg_clcf::Float64
    local_clcfs::Vector{Float64}
    wedge_counts::Vector{Int64}
    tri_counts::Vector{Int64}
end

"""
undir_clcfs
-----------------
Compute undirected closure coefficients. The input graph is given as a sparse matrix. 

A::SparseMatrixCSC{Int64,Int64}
    The adjacency matrix of an undirected graph.

returns type undir_clcfs_data
"""
function undir_clcfs(A::SparseMatrixCSC{Int64,Int64})
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    if !issymmetric(A)
        warn("Input graph is not symmetric. Symmetrizing...")
        A = max.(A, A')
    end

    n = size(A, 1)
    degs = vec(sum(A, 2))
    wedge_counts = A * (degs - 1)

    deg_order = zeros(Int64, n)
    deg_order[sortperm(degs)] = collect(1:n)
    tri_counts = zeros(Int64, n )
    for i = 1:n
        pos_i = deg_order[i]
        nbrs = [v for v in nz_row_inds(A, i) if deg_order[v] > pos_i]
        for (j, k) in combinations(nbrs, 2)
            if B[j, k] > 0; tri_counts[[i, j, k]] += 1; end
        end
    end

    local_clcfs = zeros(Float64, n)
    nz_inds = find(wedge_counts .> 0)
    local_clcfs[nz_inds] = 2 * tri_counts[nz_inds] ./ wedge_counts[nz_inds]
    return undir_clcfs_data(sum(tri_counts) / sum(degs .^ 2),
                            mean(local_clcfs[nz_inds]),
                            local_clcfs, wedge_counts, tri_counts)
end

end # module
