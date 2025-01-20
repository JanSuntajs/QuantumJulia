# ## Time evolution following ED

using DrWatson
@quickactivate :QuantumJulia

using QuantumJulia
using BenchmarkTools
using Random: rand
import PyPlot as plt

# We now select some sensible model parameters.
# Feel free to play around with them.
L::Int64= 10
nup::Int64 = 5
J::Float64 = 1.
Δ::Float64 = 1.
W::Float64 = 5.

params = @dict L nup J Δ W
fields::Vector{Float64} = W * rand(L)

basis = Basis(L, nup)
# This will create the Hamiltonian, but not build it yet:
ham = XXZHamiltonian(basis, J, Δ, fields)

@btime buildham!(ham)

println(ham._isset)

# ## Diagonalization
# We now use the full exact diagonalization to obtain
# all eigenvalues and eigenvectors of the Hamiltonian.
# Proceeed with care as this can be very resource intensive.
λ, V = diagonalize(ham; method=:dense);

# After diagonalization, we can now proceed with the time evolution,
# and, for instance, calculate the survival probability. 
# We will be interested in calculating the following quantity:
# $$
#   P(t) = |\langle \psi(0) | \psi(t) \rangle|^2
# $$
#
# Expanding, we need to evaluate the following:
# $$
#  P(t) = |\sum_\alpha |c_\alpha|^2 \exp(-i E_\alpha t)|^2.
# $$
# Each $c_\alpha$ is just the result of a dot product between the
# column of the V matrix and the initial state vector. We have
# conveniently prepared all the neccessary bits in the SurvivalProbability
# struct. To see the implementation, see `src/states.jl`.
# 
# First, we prepare the random Gaussian state:
ψ0 = random_gaussian_state(basis.nstates)
sprob = SurvivalProbability(λ, V, ψ0)
# We have now just initiated the survival probability struct but
# have not yet calculated the survival probability at a given time.
times = exp10.(range(-3, stop=4, length=1000));
probvals = sprob.(times);

# ## Saving the results using DrWatson
# See [this](https://juliadynamics.github.io/DrWatson.jl/stable/workflow/) for
# more tips and information on how to leverage DrWatson to save results.

# To format the filename, do: 
filename = savename("xxz_benchmark", params, "jld2")
# We will be using the `JLD2` format to store the data in a self-descriptive way 
# appropriate for Julia workflows. To store the results, the saving function expects 
# a dictionary of key-value pairs, with the former being the variable names and the 
# latter the results (or other quantities) to be stored.
# *Note*: if you plan to mix different programming languages, some other format, such
# as HDF5, might be more appropriate.

# This will parse the results as a dictionary which we will later
# save using the `tagsave` function.
results = @dict times probvals V λ ψ0

# DrWatson allows for many different ways to save files. This one is particularly useful,
# as it also encludes the tag of the git commit indicating the state of the project at the time
# of the simulation. If one commits the code on a regular basis, this can be very useful in
# ensuring reproducibility of the results.
tagsave(datadir("sims", filename), results);

# Of course, saving routines are of no use if the data cannot be loaded later on. Let's try
# that as well. Let's load the first simulation we saved.
firstsim = readdir(datadir("sims"))[1]
println(firstsim)
wload(datadir("sims", firstsim))
#
fig = plt.figure()

plt.loglog(times, probvals)
plt.gcf()

