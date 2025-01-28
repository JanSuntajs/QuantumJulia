# # Sparse matrices and Kronecker products
# In this notebook, we try out some sparse matrix tools for Hamiltonian
# construction based on Kronecker products. Note that this is not fully
# optimized, but it should give you an idea of how to use these tools.
# Reat [this](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) to get
# started with sparse matrices in Julia, along with the list of
# [external packages](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Noteworthy-External-Sparse-Packages).
# With the latter, make sure they are compatible with your hardware and are
# actively mantained, which is something to keep in mind when working on
# long-term projects.

using DrWatson
@quickactivate :QuantumJulia

using BenchmarkTools
using SparseArrays
import QuantumJulia.SpOperators as sop
import QuantumJulia.Operators as op
import QHam.Spin1d as sp1d
# To check the implementation of the sparse matrix approach, take a look at
# `src/operators/spoperators.jl` and `src/operators/sparseops.jl`.

# Here, let's just call the function that builds the Hamiltonian:

N = 18  # system size, no spin conservation
local_terms = [[sop.sx, sop.sy, sop.sz], [sop.sx, sop.sy, sop.sz]]  # operators at each site ()
local_couplings = [[1., 1., 1.] for _ in 0:N-2]  # couplings for each operator at each site
sitelist = [[i, (i+1)%N] for i in 0:N-2];  #  nearest neighbor interactions, two-body terms; 
#= note that the code does not handle pbc properly as some edge cases should be considered 
which is not the point of this tutorial =#
@btime H = sop.buildham(local_terms, local_couplings, sitelist, N)

# ## Comparison to a different approach
L = UInt32(18);  # system size; we keep it small here, but do play around with it
nup = nothing;  # number of up spins for a particle-number conserving model

#= True for whether the state indices are stored or not. In the opposite case, binary search
   is performed to find the index of the state in the full basis.=#
basis = sp1d.SpinBasis(L, nup, false)
coups = [[1., i, i+1] for i in 0:L-2];

input_list = [["xx", coups], ["yy", coups], ["zz", coups]]

ham = sp1d.Ham(basis, input_list, ComplexF64);

sp1d.buildham!(ham);

# ## Some benchmarks
@benchmark H=sop.buildham(local_terms, local_couplings, sitelist, N)
#
@benchmark sp1d.buildham!(ham) setup=(ham = sp1d.Ham(basis, input_list, ComplexF64))

# ### Out of curiosity
# Let's see the construction if lookup table is used instead of binary search
basis = sp1d.SpinBasis(L, nup, true)

ham = sp1d.Ham(basis, input_list, ComplexF64);

@benchmark sp1d.buildham!(ham) setup=(ham = sp1d.Ham(basis, input_list, ComplexF64))

# Bottom line: Sparse matrix approac via Kronecker products is a worthwile avenue to explore,
# but I'd need to get into the details to optimize it even further. Also, there's a bit of a
# conundrum with the high-performance sparse options.