module ClosureCoefficients

using Combinatorics

export undir_clcfs, undir_clcfs_data

const WEDGE_OO = 1
const WEDGE_OI = 2
const WEDGE_IO = 3
const WEDGE_II = 4

const TRI_OO_O = 1
const TRI_OO_I = 2
const TRI_OI_O = 3
const TRI_OI_I = 4
const TRI_IO_O = 5
const TRI_IO_I = 6
const TRI_II_O = 7
const TRI_II_I = 8

nz_row_inds(A::SparseMatrixCSC{Int64,Int64}, ind::Int64) =
    A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

"""
clcfs_data
----------
Closure coefficients data.

The field values are
global_clcf::Float64
    The global closure coefficient.
avg_clcf::Float64 
    The average closure coefficient (the mean is taken over nodes that head at
    least one wedge).
local_clcfs::Vector{Float64}
    Vctor of local closure coefficients. If a node is not the head of at least
    one wedge, then the value is 0.
wedges::Vector{Int64}
    Vector of the number of wedges at each node: wedges[v] is the number of
    wedges headed by node v.
triangles::Vector{Int64}
    Vector of triangle counts of nodes: tri_counts[v] is the number of triangles
    containing v.

"""
immutable clcfs_data
    global_clcf::Float64
    avg_clcf::Float64
    local_clcfs::Vector{Float64}
    wedges::Vector{Int64}
    triangles::Vector{Int64}
end

"""
dir_clcfs
-----------------
Compute directed closure coefficients. The input graph is given as a sparse
matrix A, where A[i, j] = 1 if i --> j in the graph.

A::SparseMatrixCSC{Int64,Int64}
    The adjacency matrix of an undirected graph.

returns type clcfs_data
"""
function dir_clcfs(A::SparseMatrixCSC{Int64,Int64})
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    At = A'
    B = min.(A, A')

    n = size(A, 1)
    dB = vec(sum(B, 1))
    dI = vec(sum(A, 1))
    dO = vec(sum(At, 1))
    
    wedges_OO = A  * dO - dB
    wedges_IO = At * (dO - 1)
    wedges_OI = A  * (dI - 1)
    wedges_II = At * dI - dB

    C = max.(A, A')
    n = size(C, 1)
    deg_order = zeros(Int64, n)
    deg_order[sortperm(vec(sum(C, 1)))] = collect(1:n)
    triangles = zeros(Int64, 8, n)
    for i = 1:n
        pos_i = deg_order[i]
        nbrs = [v for v in nz_row_inds(C, i) if deg_order[v] > pos_i]
        for (j, k) in combinations(nbrs, 2)
            if C[j, k] > 0
                # Triangle between i, j, k
                Aij = A[i, j] > 0
                Aik = A[i, k] > 0
                Aji = A[j, i] > 0
                Ajk = A[j, k] > 0
                Aki = A[k, i] > 0
                Akj = A[k, j] > 0
                if Aij && Ajk && Aik
                    triangles[TRI_OO_O, i] += 1
                    triangles[TRI_OI_O, i] += 1
                    triangles[TRI_OI_I, j] += 1
                    triangles[TRI_IO_O, j] += 1
                    triangles[TRI_IO_I, k] += 1
                    triangles[TRI_II_I, k] += 1                    
                end
                if Aij && Ajk && Aki
                    triangles[TRI_OO_I, i] += 1
                    triangles[TRI_II_O, i] += 1
                    triangles[TRI_OO_I, j] += 1
                    triangles[TRI_II_O, j] += 1
                    triangles[TRI_OO_I, k] += 1
                    triangles[TRI_II_O, k] += 1                    
                end
                if Aij && Akj && Aik
                    triangles[TRI_OI_O, i] += 1
                    triangles[TRI_OO_O, i] += 1
                    triangles[TRI_II_I, j] += 1
                    triangles[TRI_IO_I, j] += 1
                    triangles[TRI_IO_O, k] += 1
                    triangles[TRI_OI_I, k] += 1                    
                end
                if Aij && Akj && Aki
                    triangles[TRI_OI_I, i] += 1
                    triangles[TRI_IO_O, i] += 1
                    triangles[TRI_IO_I, j] += 1
                    triangles[TRI_II_I, j] += 1
                    triangles[TRI_OO_O, k] += 1
                    triangles[TRI_OI_O, k] += 1                    
                end
                if Aji && Ajk && Aik
                    triangles[TRI_IO_O, i] += 1
                    triangles[TRI_OI_I, i] += 1
                    triangles[TRI_OI_O, j] += 1
                    triangles[TRI_OO_O, j] += 1
                    triangles[TRI_II_I, k] += 1
                    triangles[TRI_IO_I, k] += 1                    
                end
                if Aji && Ajk && Aki
                    triangles[TRI_IO_I, i] += 1
                    triangles[TRI_II_I, i] += 1
                    triangles[TRI_OO_O, j] += 1
                    triangles[TRI_OI_O, j] += 1
                    triangles[TRI_OI_I, k] += 1
                    triangles[TRI_IO_O, k] += 1                    
                end
                if Aji && Akj && Aik
                    triangles[TRI_II_O, i] += 1
                    triangles[TRI_OO_I, i] += 1
                    triangles[TRI_II_O, j] += 1
                    triangles[TRI_OO_I, j] += 1
                    triangles[TRI_II_O, k] += 1
                    triangles[TRI_OO_I, k] += 1                    
                end
                if Aji && Akj && Aki
                    triangles[TRI_II_I, i] += 1
                    triangles[TRI_IO_I, i] += 1
                    triangles[TRI_IO_O, j] += 1
                    triangles[TRI_OI_I, j] += 1
                    triangles[TRI_OI_O, k] += 1
                    triangles[TRI_OO_O, k] += 1                    
                end
            end
        end
    end
    triangles = triangles'

    function create_data(triangles::Vector{Int64}, wedges::Vector{Int64})
        local_clcfs = zeros(Float64, n)
        nz_inds = find(wedges .> 0)
        local_clcfs[nz_inds] = triangles[nz_inds] ./ wedges[nz_inds]
        return clcfs_data(sum(tri_counts) / sum(wedges),
                          mean(local_clcfs[nz_inds]), local_clcfs, wedges, triangles)
    end
    ret = Vector{clcfs_data}(8)
    ret[TRI_OO_O] = create_data(triangles[:, TRI_OO_O], wedges_OO)
    ret[TRI_OO_I] = create_data(triangles[:, TRI_OO_I], wedges_OO)    
    ret[TRI_OI_O] = create_data(triangles[:, TRI_OI_O], wedges_OI)
    ret[TRI_OI_I] = create_data(triangles[:, TRI_OI_I], wedges_OI)    
    ret[TRI_IO_O] = create_data(triangles[:, TRI_IO_O], wedges_IO)
    ret[TRI_IO_I] = create_data(triangles[:, TRI_IO_I], wedges_IO)    
    ret[TRI_II_O] = create_data(triangles[:, TRI_II_O], wedges_II)
    ret[TRI_II_I] = create_data(triangles[:, TRI_II_I], wedges_II)    
    return ret
end


"""
undir_clcfs
-----------------
Compute undirected closure coefficients. The input graph is given as a symmetric
sparse matrix.

A::SparseMatrixCSC{Int64,Int64}
    The adjacency matrix of an undirected graph.

returns type clcfs_data
"""
function undir_clcfs(A::SparseMatrixCSC{Int64,Int64})
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    if !issymmetric(A)
        warn("Input graph is not symmetric. Symmetrizing...")
        A = max.(A, A')
    end

    n = size(A, 1)
    degs = vec(sum(A, 1))
    wedges = A * (degs - 1)

    deg_order = zeros(Int64, n)
    deg_order[sortperm(degs)] = collect(1:n)
    triangles = zeros(Int64, n)
    for i = 1:n
        pos_i = deg_order[i]
        nbrs = [v for v in nz_row_inds(A, i) if deg_order[v] > pos_i]
        for (j, k) in combinations(nbrs, 2)
            if A[j, k] > 0; triangles[[i, j, k]] += 1; end
        end
    end

    local_clcfs = zeros(Float64, n)
    nz_inds = find(wedges .> 0)
    local_clcfs[nz_inds] = 2 * triangles[nz_inds] ./ wedges[nz_inds]
    return clcfs_data(2 * sum(triangles) / sum(wedges),
                      mean(local_clcfs[nz_inds]), local_clcfs, wedges, triangles)
end

end # module
