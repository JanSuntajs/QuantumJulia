"""
In this script, we benchmark the performance of the `buildham!` function.
To run as a script, do the following:

```bash
julia --project=. scripts/01_simple_xxz.jl
```
Here, we have assumed you are running the code from the root folder of the
project. Since we are managing the project using DrWatson, this is not strictly
needed. Even running from elsewhere (and without the --project=.) line will be fine
as the `@quickactivate` macro will take care of the environment activation.

"""

using DrWatson
@quickactivate :QuantumJulia

using QuantumJulia
using BenchmarkTools
using Random: rand

# Feel free to update and ramp these parameters up
const L::Int64 = 16
const nup::Int64 = 8
const J::Float64 = 1.
const Δ::Float64 = 1.
const W::Float64 = 5.
const fields::Vector{Float64} = W * rand(L)

params = @dict L nup J Δ W

const basis::Basis = Basis(L, nup)
# This will create the Hamiltonian, but not build it yet.
ham::XXZHamiltonian = XXZHamiltonian(basis, J, Δ, fields)

@btime buildham!(ham)

println(ham._isset)

# Diagonalize
# we use the method eigsolve from KrylovKit in this
# case. We use the default settings of the solver, but the
# formulation of diagonalize allows for both the specification of
# the additional positional and keyword arguments
# https://jutho.github.io/KrylovKit.jl/stable/man/eig/
println("Starting diagonalization!")
@btime λ, V = diagonalize(ham; method=:dense);
println("Diagonalization done!")
λ, V = diagonalize(ham; method=:dense);
# Now calculate the survival probability.
ψ0 = random_gaussian_state(basis.nstates)
@time sprob = SurvivalProbability(λ, V, ψ0)

times = exp10.(range(-3, stop=4, length=1000));
probvals = similar(times)
@time probvals .= sprob.(times);

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
