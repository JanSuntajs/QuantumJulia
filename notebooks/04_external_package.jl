# # Using external packages
#
# In this notebook, we take a look at using external packages in our Julia code.
# Specifically, we will make use of the `QHam.jl` package which I developed myself
# to handle the Hamiltonian construction and diagonalization in quantum many-body
# systems as well as non-interacting models. The package is not yet publicly
# registedred, so this is a good chance to take a look at how to use unregistered
# packages in Julia. Note that the tutorial is to a large extent based
# on the tutorials for the package.
#
# ## Installation
#
# ```julia
# #  in Pkg mode entered by pressing ]
# add git@github.com:JanSuntajs/QHam.jl.git  # note you need appripriate permissions
# ```
#
# ## Usage
# Below, we show how to use the package in the context of the simple XXZ model.
using DrWatson
@quickactivate :QuantumJulia

using Random
using BenchmarkTools
import PyPlot as plt

import QHam.Hamiltonian: Ham
import QHam.Diagonalize: diagonalize, diagonalize!
using QHam.Models.Spin1d

# ## Building the parameter list and defining the Hamiltonian

# Here, we'll show how to prepare 1D hamiltonian, for instance that of the [XXZ model](https://arxiv.org/pdf/2004.01719.pdf).
# The Hamiltonian reads:

# $$
# \hat{H} = J \sum\limits_{l=1}^L \left(\hat{s}^x_{l}\hat{s}^x_{l+1}+ \hat{s}^y_{l}\hat{s}^y_{l+1} + \Delta_l\hat{s}^z_{l}\hat{s}^z_{l+1}\right) + \sum\limits_{l=1}^L w_l \hat{s}^z_l,
# $$
# where $\hat{s}^\alpha_l$ are the spin-1/2 operators at site $l$ and $\alpha=(x, y, z)$ denotes the spin 
# operator type. To define the Hamiltonian in our program, we'll first define the model constants, as shown below.

L = UInt32(12);  # system size; we keep it small here, but do play around with it
nup = UInt32(L // 2);  # number of up spins for a particle-number conserving model
J = 1.0;  # Heisenberg exchange
Δ = 1.0;  # anisotropy parameter

# ## Defining the basis
# Next, we can define the Hamiltonian basis which is a neccesary step in the Hamiltonian 
# construction - we first need to construct the basis on which the Hamiltonian acts to 
# obtain its matrix representation. Below, we construct the basis for a given $L$ and `nup` spins. Note that 
# if `nup = nothing`, the spin basis for the full Hilbert space is constructed.
# Apart from constructing the basis, we elso show how to benchmark the process of basis construction for performance
# using julia built-in tools.

#= True for whether the state indices are stored or not. In the opposite case, binary search
   is performed to find the index of the state in the full basis.=#
basis = SpinBasis(L, nup, true)
# We can also take a look at the documentation of the `SpinBasis` struct:
@doc SpinBasis(L, nothing)
# To see the number of basis states, we do:
Int(basis.nstates)

# ## Explaining the parameter lists
# We now prepare the parameter lists to build the Hamiltonian. The logic behind our approach is very 
# close to the one in the Python package [quspin](http://quspin.github.io/QuSpin/). Namely, 
# for each term in the Hamiltonian, we define a *coupling list*, which defines couplings between different sites
# or fields acting on individual sites. We explain this by a commented example below:

# below, we define the couplings in the XXZ Hamiltonian. The structure of the
# terms is standard and is as follows:
# we provide a vector of vectors, which we'll denote as a coupling vector.
# In each of its inner vectors the structure is standard, with the first value
# denoting the 
# value of the coupling/field. The remainder of the inner vector denotes the site(s)
# on which the operator term acts. Below, we prepare the terms for the 
# +- term, the interaction term and the random fields. Note that in the latter, 
# we only provide one value for the site, since the fields are represented by
# single body operators.
# *NOTE ALSO* #1: we start our indexing with zero and use the i%L syntax to implement
# the periodic boundary conditions
# *NOTE ALSO* #2: we explain how to couple this information with the actual operators below =#
spsm = [[J, i, (i + 1) % L] for i in 0:L-1];
inter = [[J * Δ, i, (i + 1) % L] for i in 0:L-1];
fields = [[(0.1 / L) * (-1.0)^i, i] for i in 0:L-1];

# To define a single Hamiltonian term, we pass a mixed-type vector:
# its first entry denotes the operator string and its second entry is 
# the coupling vector. 

# In passing the operator string, one needs to take care that it's entries
# are among the allowed operators for a given model type. Also, the length
# of the operator string should match the number of sites appearing in the
# inner vectors of the coupling vectors. In other words, since "+-" defines
# a two-site operator, the corresponding coupling vector should be defined as
# [[J, i, (i + 1) % L] for i in 0:L-1]. If not, the program will throw an error.

# Below, we compose the "+-" and "-+" operators, as well as the interaction term
# "zz" and the fields "z". Note that there is technically no limit to the number
# of sites one could couple together in composing a N-body operator. 

input_list = [["zz", inter], ["+-", spsm], ["-+", spsm], ["z", fields]]

# ## Defining the Hamiltonian and building its matrix
ham = Ham(basis, input_list, Float64);
#= check if the Hamiltonian matrix is set or not =#
println(ham.isset)
buildham!(ham);
#= check for the matrix again =#
println(ham.isset)
#= print the Ham docs =#
@doc ham

# Benchmark the process of building the Hamiltonian matrix
@btime buildham!(ham)

# ## Visualization
# To visualize the matrix, we need to access the `ham.mat` field.

plotmat = Array(abs.(ham.mat))
cm = plt.get_cmap(:gray_r)

fig, ax = plt.subplots()

h = plt.imshow(plotmat, cm, vmin=0.0, vmax=maximum(plotmat))
plt.colorbar(h, ax=ax)
display(fig)

# ## Diagonalization
# Already at the package level (and without the fancy shift and invert or Polfed
# implementations), we can diagonalize the Hamiltonian using the `diagonalize` or
# `diagonalize!` routines. 

# To perform full exact diagonalization and return the eigenvalues and eigenvectors,
# we do the following:

λ, V = diagonalize(ham; method=:dense);

# To perform the same operation in place (after which the eigenvalues and eigenvectors
# are stored in the `ham2` object), we do:

diagonalize!(ham; method=:dense)

# We can also use the KrylovKit package to perform the diagonalization. For that, we
# use the `:krylov` method. The `eigsolve` function from the KrylovKit package is used

λ, V = diagonalize(ham, 5, :LM; method=:krylov)

# Let us also take a look at the documentation of the diagonalization routines for some
# additional information:

@doc diagonalize

#
@doc diagonalize!