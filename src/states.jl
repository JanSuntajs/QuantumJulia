"""
    random_gaussian_state(rng::AbstractRNG, dim::Integer) -> Vector{ComplexF64}

Generates a random Gaussian quantum state with the specified dimension.

# Arguments
- `rng::AbstractRNG`: The random number generator to use.
- `dim::Integer`: The dimension of the quantum state.

# Returns
- `Vector{ComplexF64}`: A normalized random Gaussian quantum state.

# Details
The coefficients of the quantum state are drawn from a complex Gaussian 
distribution with zero mean and unit variance. The state is then normalized 
to ensure it represents a valid quantum state.

# Example
```julia
rng = MersenneTwister(1234)
ψ = random_gaussian_state(rng, 10)
println("Random Gaussian state: ", ψ)
println("Norm of the state: ", norm(ψ))
```
"""
function random_gaussian_state(dim::T, 
    ::Type{M}=ComplexF64,
    rng::AbstractRNG=MersenneTwister()) where {T <: Integer, M <: Complex}
    ψ = randn(rng, M, dim)
    ψ ./= norm(ψ)
    return ψ
end

"""
    struct SurvivalProbability{T <: Real}

A struct representing the survival probability of a quantum state over time.

# Fields
- `λ::Vector{T}`: The eigenvalues of the Hamiltonian.
- `coeffs::Vector{Complex{T}}`: The coefficients of the initial state in the eigenbasis.

# Constructor
- `SurvivalProbability(λ::Vector{T}, U::AbstractMatrix{Complex{T}}, ψ0::AbstractVector{Complex{T}}; coeffs=similar(λ, Complex{T}))`: Creates a new SurvivalProbability instance.

# Example
```julia
λ = rand(10)
U = rand(ComplexF64, 10, 10)
ψ0 = rand(ComplexF64, 10)
sprob = SurvivalProbability(λ, U, ψ0)
```
"""
struct SurvivalProbability{T <: Real}
    λ::Vector{T}
    coeffs::Vector{Complex{T}}

    function SurvivalProbability(λ::Vector{T}, 
        U::AbstractMatrix{M}, 
        ψ0::AbstractVector{Complex{T}}) where {T <: Real, M <: Union{T, Complex{T}}}

        coeffs = abs2.(U' * ψ0)

        new{T}(λ, coeffs)
    end
end

"""
    (sprob::SurvivalProbability)(t::Real) -> Complex{T}

Evaluates the survival probability at time `t`.

# Arguments
- `sprob::SurvivalProbability`: The SurvivalProbability instance.
- `t::Real`: The time at which to evaluate the survival probability.

# Returns
- `Complex{T}`: The survival probability at time `t`.

# Example
```julia
sprob = SurvivalProbability(rand(10), rand(ComplexF64, 10, 10), rand(ComplexF64, 10))
prob = sprob(1.0)
```
"""
@inline function (sprob::SurvivalProbability{T})(t::Real)::T  where {T <: Real}
    return abs2(sprob.coeffs ⋅ exp.(-1im * sprob.λ * t))
end

