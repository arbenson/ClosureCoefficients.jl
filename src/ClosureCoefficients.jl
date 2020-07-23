module ClosureCoefficients

using Combinatorics
using LinearAlgebra
using Printf
using SparseArrays
using StatsBase

export undir_clcfs, dir_clcfs, clcfs_data, construct_directed_data, load_example_data

const SpIntMat = SparseMatrixCSC{Int64,Int64}

const TRI_oo_o = 1
const TRI_oo_i = 2
const TRI_oi_o = 3
const TRI_oi_i = 4
const TRI_io_o = 5
const TRI_io_i = 6
const TRI_ii_o = 7
const TRI_ii_i = 8

"""
Read a graph from a text file. Skips any lines that begin with '#' or '%'
characters. Nodes are mapped to values 1, ..., n.

Returns (A, v), where
   A is the n x n undirected graph adjacency matrix
   v is a length-n vector such that v[i] is the original node identifier in the text file
"""
function read_graph_txt(filename::AbstractString)
    # index mapping
    index_map = Dict{Int64,Int64}()
    index_map_vec = Int64[]
    function get_mapped_index(x::Int64)
        if !haskey(index_map, x)
            next_index = length(index_map) + 1
            index_map[x] = next_index
            push!(index_map_vec, x)
            return next_index
        end
        return index_map[x]
    end

    # Read data
    I = Int64[]
    J = Int64[]
    open(filename) do f
        for line in eachline(f)
            # Skip lines starting with '#' or '%'
            if line[1] == '#' || line[1] == '%'; continue; end
            edge = split(line)
            u = parse(Int64, edge[1])
            v = parse(Int64, edge[2])
            push!(I, get_mapped_index(u))
            push!(J, get_mapped_index(v))
        end
    end

    n = max(maximum(I), maximum(J))
    A = convert(SpIntMat, sparse(I, J, ones(length(I)), n, n))    
    return (A, index_map_vec)
end

"""
Load an example file from the data directory.
"""
function load_example_data(filename::AbstractString; symm::Bool=false)
    pathname = joinpath(dirname(dirname(@__FILE__)), "data")
    filename = joinpath(pathname, filename)
    if   isfile(filename)
        A = read_graph_txt(filename)[1]
        if symm
            A = max.(A, A')
        end
        return A
    else
        error(@sprintf("Could not find file %s", name))
    end
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
struct clcfs_data
    global_clcf::Float64
    avg_clcf::Float64
    local_clcfs::Vector{Float64}
    wedges::Vector{Int64}
    triangles::Vector{Int64}
end

function construct_directed_data(triangles::Vector{Int64}, wedges::Vector{Int64})
    local_clcfs = zeros(Float64, length(wedges))
    nz_inds = findall(wedges .> 0)
    local_clcfs[nz_inds] = triangles[nz_inds] ./ wedges[nz_inds]
    avg_clcf = (length(nz_inds) > 0) ? mean(local_clcfs[nz_inds]) : 0.0
    num_wedges = sum(wedges)
    global_clcf = (num_wedges > 0) ? sum(triangles) / num_wedges : 0.0
    return clcfs_data(global_clcf, avg_clcf, local_clcfs, wedges, triangles)
end

function degree_order(A::SpIntMat)
    n = size(A, 1)
    deg_order = zeros(Int64, n)
    deg_order[sortperm(vec(sum(A, dims=1)))] = collect(1:n)
end

nz_row_inds(A::SpIntMat, ind::Int64) =
    A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

# Return generator for pairs of neighbors higher in the ordering
function ordered_nbr_pairs(A::SpIntMat, i::Int64, order::Vector{Int64})
    pos_i = order[i]
    return combinations(filter(nbr -> order[nbr] > pos_i, nz_row_inds(A, i)), 2)
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
    A -= Diagonal(A)
    At = convert(SpIntMat, A')
    B = min.(A, At)

    n = size(A, 1)
    dB = vec(sum(B, dims=1))
    dI = vec(sum(A, dims=1))
    dO = vec(sum(At, dims=1))

    @show typeof(A), typeof(dO), typeof(A * dO)
    
    wedges_oo = A  * dO - dB
    wedges_io = At * (dO .- 1)
    wedges_oi = A  * (dI .- 1)
    wedges_ii = At * dI - dB

    C = max.(A, At)
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
                    triangles[TRI_oo_o, i] += 1
                    triangles[TRI_oi_o, i] += 1
                    triangles[TRI_oi_i, j] += 1
                    triangles[TRI_io_o, j] += 1
                    triangles[TRI_io_i, k] += 1
                    triangles[TRI_ii_i, k] += 1
                end
                if Aij && Ajk && Aki
                    triangles[TRI_oo_i, i] += 1
                    triangles[TRI_ii_o, i] += 1
                    triangles[TRI_oo_i, j] += 1
                    triangles[TRI_ii_o, j] += 1
                    triangles[TRI_oo_i, k] += 1
                    triangles[TRI_ii_o, k] += 1
                end
                if Aij && Akj && Aik
                    triangles[TRI_oi_o, i] += 1
                    triangles[TRI_oo_o, i] += 1
                    triangles[TRI_ii_i, j] += 1
                    triangles[TRI_io_i, j] += 1
                    triangles[TRI_io_o, k] += 1
                    triangles[TRI_oi_i, k] += 1
                end
                if Aij && Akj && Aki
                    triangles[TRI_oi_i, i] += 1
                    triangles[TRI_io_o, i] += 1
                    triangles[TRI_io_i, j] += 1
                    triangles[TRI_ii_i, j] += 1
                    triangles[TRI_oo_o, k] += 1
                    triangles[TRI_oi_o, k] += 1
                end
                if Aji && Ajk && Aik
                    triangles[TRI_io_o, i] += 1
                    triangles[TRI_oi_i, i] += 1
                    triangles[TRI_oi_o, j] += 1
                    triangles[TRI_oo_o, j] += 1
                    triangles[TRI_ii_i, k] += 1
                    triangles[TRI_io_i, k] += 1
                end
                if Aji && Ajk && Aki
                    triangles[TRI_io_i, i] += 1
                    triangles[TRI_ii_i, i] += 1
                    triangles[TRI_oo_o, j] += 1
                    triangles[TRI_oi_o, j] += 1
                    triangles[TRI_oi_i, k] += 1
                    triangles[TRI_io_o, k] += 1
                end
                if Aji && Akj && Aik
                    triangles[TRI_ii_o, i] += 1
                    triangles[TRI_oo_i, i] += 1
                    triangles[TRI_ii_o, j] += 1
                    triangles[TRI_oo_i, j] += 1
                    triangles[TRI_ii_o, k] += 1
                    triangles[TRI_oo_i, k] += 1
                end
                if Aji && Akj && Aki
                    triangles[TRI_ii_i, i] += 1
                    triangles[TRI_io_i, i] += 1
                    triangles[TRI_io_o, j] += 1
                    triangles[TRI_oi_i, j] += 1
                    triangles[TRI_oi_o, k] += 1
                    triangles[TRI_oo_o, k] += 1
                end
            end
        end
    end
    triangles = triangles'

    ret = Dict{String, clcfs_data}()
    ret["oo_o"] = construct_directed_data(triangles[:, TRI_oo_o], wedges_oo)
    ret["oo_i"] = construct_directed_data(triangles[:, TRI_oo_i], wedges_oo)
    ret["oi_o"] = construct_directed_data(triangles[:, TRI_oi_o], wedges_oi)
    ret["oi_i"] = construct_directed_data(triangles[:, TRI_oi_i], wedges_oi)
    ret["io_o"] = construct_directed_data(triangles[:, TRI_io_o], wedges_io)
    ret["io_i"] = construct_directed_data(triangles[:, TRI_io_i], wedges_io)
    ret["ii_o"] = construct_directed_data(triangles[:, TRI_ii_o], wedges_ii)
    ret["ii_i"] = construct_directed_data(triangles[:, TRI_ii_i], wedges_ii)
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
    A -= Diagonal(A)
    if !issymmetric(A)
        println("Input graph is not symmetric. Symmetrizing...")
        A = max.(A, A')
    end

    n = size(A, 1)
    degs = vec(sum(A, dims=1))
    wedges = A * (degs .- 1)

    deg_order = degree_order(A)
    triangles = zeros(Int64, n)
    for i = 1:n
        for (j, k) in ordered_nbr_pairs(A, i, deg_order)
            if A[j, k] > 0; triangles[[i, j, k]] .+= 1; end
        end
    end

    local_clcfs = zeros(Float64, n)
    nz_inds = findall(wedges .> 0)
    local_clcfs[nz_inds] = 2 * triangles[nz_inds] ./ wedges[nz_inds]
    return clcfs_data(2 * sum(triangles) / sum(wedges),
                      mean(local_clcfs[nz_inds]), local_clcfs, wedges, triangles)
end

end # module
