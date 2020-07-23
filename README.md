# ClosureCoefficients.jl

This package provides an interface to compute "closure coefficients" in networks, an idea developed in the following papers:

- [Measuring Directed Triadic Closure with Closure Coefficients](http://www.cs.cornell.edu/~arb/papers/directed-closure-NWS-2020.pdf). Hao Yin, Austin R. Benson, and Johan Ugander. Network Science, 2020.
- [The Local Closure Coefficient: A New Perspective On Network Clustering](http://www.cs.cornell.edu/~arb/papers/closure-coefficients-WSDM-2019.pdf). Hao Yin, Austin R. Benson, and Jure Leskovec. Proceedings of the ACM International Conference on Web Search and Data Mining (WSDM), 2019.

## Installation
From Julia
```julia
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/arbenson/ClosureCoefficients.jl"))
```

## Examples
The following examples assume that you have added the package and run the following in Julia.
```julia
using ClosureCoefficients
```

#### Compute directed clousure coefficients of the Florida Bay food web.
```julia
A = load_example_data("FW-Florida.txt", symm=false)
clcfs = dir_clcfs(A)
clcfs["ii_o"].avg_clcf  # mean directed closure for o-closed ii wedges
clcfs["oi_i"].local_clcfs  # fraction i-closed oi wedges, per node
```

A "zero" value for a local closure coefficient can mean two things: the node is at the head of at least one wedge but the wedge never closes or (ii) the node is not the head of a wedge. It is easy to find the nodes in the latter case because the data structure also returns the wedge counts.
```julia
findall(clcfs["oi_i"].wedges .== 0)
```

#### Compute undirected clousure coefficients on the arXiv-AstroPh network.

```julia
A = load_example_data("arxiv-AstroPh.txt", symm=true)
clcfs = undir_clcfs(A)
clcfs.avg_clcf
```