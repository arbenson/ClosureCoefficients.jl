module ClosureCoefficients

using Combinatorics

export undir_clcfs, dir_clcfs, clcfs_data

const SpIntMat = SparseMatrixCSC{Int64,Int64}

const TRI_OO_O = 1
const TRI_OO_I = 2
const TRI_OI_O = 3
const TRI_OI_I = 4
const TRI_IO_O = 5
const TRI_IO_I = 6
const TRI_II_O = 7
const TRI_II_I = 8

nz_row_inds(A::SpIntMat, ind::Int64) =
    A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

function create_data_dir(triangles::Vector{Int64}, wedges::Vector{Int64})
    local_clcfs = zeros(Float64, length(wedges))
    nz_inds = find(wedges .> 0)
    local_clcfs[nz_inds] = triangles[nz_inds] ./ wedges[nz_inds]
    avg_clcf = (length(nz_inds) > 0) ? mean(local_clcfs[nz_inds]) : 0.0
    num_wedges = sum(wedges)
    global_clcf = (num_wedges > 0) ? sum(triangles) / num_wedges : 0.0
    return clcfs_data(global_clcf, avg_clcf, local_clcfs, wedges, triangles)
end

function degree_order(A::SpIntMat)
    n = size(A, 1)
    deg_order = zeros(Int64, n)
    deg_order[sortperm(vec(sum(A, 1)))] = collect(1:n)
end

# Return generator for pairs of neighbors higher in the ordering
function ordered_nbr_pairs(A::SpIntMat, i::Int64, order::Vector{Int64})
    pos_i = order[i]
    return combinations(filter(nbr -> order[nbr] > pos_i, nz_row_inds(A, i)), 2)
end

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
function dir_clcfs(A::SpIntMat)
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
    deg_order = degree_order(C)
    triangles = zeros(Int64, 8, n)
    for i = 1:n
        for (j, k) in ordered_nbr_pairs(C, i, deg_order)
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

    ret = Dict{String, clcfs_data}()
    ret["OO_O"] = create_data_dir(triangles[:, TRI_OO_O], wedges_OO)
    ret["OO_I"] = create_data_dir(triangles[:, TRI_OO_I], wedges_OO)
    ret["OI_O"] = create_data_dir(triangles[:, TRI_OI_O], wedges_OI)
    ret["OI_I"] = create_data_dir(triangles[:, TRI_OI_I], wedges_OI)
    ret["IO_O"] = create_data_dir(triangles[:, TRI_IO_O], wedges_IO)
    ret["IO_I"] = create_data_dir(triangles[:, TRI_IO_I], wedges_IO)
    ret["II_O"] = create_data_dir(triangles[:, TRI_II_O], wedges_II)
    ret["II_I"] = create_data_dir(triangles[:, TRI_II_I], wedges_II)
    return ret
end

"""
dir_clcfs2
-----------------
Compute directed closure coefficients. The input graph is given as a sparse
matrix A, where A[i, j] = 1 if i --> j in the graph.

A::SparseMatrixCSC{Int64,Int64}
    The adjacency matrix of an undirected graph.

returns type clcfs_data
"""
function dir_clcfs2(A::SpIntMat)
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
    deg_order = degree_order(C)
    triangles = zeros(Int64, 8, n)
    for i = 1:n
        for (j, k) in ordered_nbr_pairs(C, i, deg_order)
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

    ret = Dict{String, clcfs_data}()
    ret["OO_O"] = create_data_dir(triangles[:, TRI_OO_O], wedges_OO)
    ret["OO_I"] = create_data_dir(triangles[:, TRI_OO_I], wedges_OO)
    ret["OI_O"] = create_data_dir(triangles[:, TRI_OI_O], wedges_OI)
    ret["OI_I"] = create_data_dir(triangles[:, TRI_OI_I], wedges_OI)
    ret["IO_O"] = create_data_dir(triangles[:, TRI_IO_O], wedges_IO)
    ret["IO_I"] = create_data_dir(triangles[:, TRI_IO_I], wedges_IO)
    ret["II_O"] = create_data_dir(triangles[:, TRI_II_O], wedges_II)
    ret["II_I"] = create_data_dir(triangles[:, TRI_II_I], wedges_II)
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
function undir_clcfs(A::SpIntMat)
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    if !issymmetric(A)
        warn("Input graph is not symmetric. Symmetrizing...")
        A = max.(A, A')
    end

    n = size(A, 1)
    degs = vec(sum(A, 1))
    wedges = A * (degs - 1)

    deg_order = degree_order(A)
    triangles = zeros(Int64, n)
    for i = 1:n
        for (j, k) in ordered_nbr_pairs(A, i, deg_order)
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
