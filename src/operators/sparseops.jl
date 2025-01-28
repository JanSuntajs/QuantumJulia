
# using BenchmarkTools

const sx = ls.staticize(sparse([1, 2], [2, 1], [0.5, 0.5], 2, 2))
const sy = ls.staticize(sparse([1, 2], [2, 1], [-0.5im, 0.5im], 2, 2))
const sz = ls.staticize(sparse([1, 2], [1, 2], [0.5, -0.5], 2, 2))
const sp = ls.staticize(sparse([2], [1], [1.0], 2, 2))  # S+ operator
const sm = ls.staticize(sparse([1], [2], [1.0], 2, 2) ) # S- operator
const id = ls.staticize(sparse(1. * I(2)))
const id1 = ls.staticize(sparse(I(1)))

const ls_CSC = ls.SSparseMatrixCSC


    # embed_operators(opvec::Vector, sitevec::Vector{Int}, N::Int)
"""
    embed_operators(opvec::Vector, sitevec::Vector{Int}, N::Int)

Embed a vector of operators into a larger Hilbert space.

# Arguments
- `opvec::Vector`: A vector of operators to embed.
- `sitevec::Vector{Int}`: A vector of sites where the operators are to be embedded.
- `N::Int`: The total number of sites in the Hilbert space.

# Returns
- A sparse matrix representing the embedded operators in the larger Hilbert space.

This function takes a vector of operators and embeds them into a larger Hilbert space of size `N`. 
The operators are placed at the specified sites in `sitevec`, with identities filling the remaining spaces.
"""
function embed_operators(opvec::Vector{T}, sitevec::Vector{Int}, N::Int) where {T <: ls_CSC}
    # Expand sitevec with leading 1 and trailing N + 1
    diffs = _make_diffs(sitevec, N)
    # Define the starting term as identity of 1
    temp = 1.#sparse(I(1))

    # Iteratively apply operators and identities
    # each operator is applied to the previous operator
    # pairs of identities and trailing operators are treated
    # together
    @inbounds for (i, op) in enumerate(opvec)
        temp_ = ls.kron(ls.IMatrix(diffs[i]), op )  
        temp = ls.kron(temp, temp_) 
        # println(size(ret))                  # Add the operator
        # ret = kron(ret, sparse(I(2^(diffs[i + 1]))))  # Add trailing identity
    end
    # the final trailing identity is treated separately
    # if sitevec[end] != N
    temp = ls.kron(temp, ls.IMatrix(diffs[end]))  # Add trailing identity
    # end

    return temp
end




function _build_local_term(local_terms, local_couplings::Vector{LC}, sitevec::Vector{SV}, N::SV) where {LC <: Number, SV <:Integer}

    diffs = _make_diffs(sitevec, N)
    temp = Ref(id1)
    # local_terms[1] .=

    @inbounds for (i, term) in enumerate(local_terms)
        temp_ = ls.kron.(Ref(ls.IMatrix(diffs[i])), term)
        temp = ls.kron.(temp, temp_)
    end
    
    temp = sum(local_couplings .* ls.kron.(temp, Ref(ls.IMatrix(diffs[end]))))
    return temp
end


function buildham(local_terms, local_couplings::Vector{Vector{LC}}, sites::Vector{Vector{SV}}, N::SV) where {LC <: Number, SV <:Integer}

    mat = ls.spzeros(2^N, 2^N)

    for (i, sitevec) in enumerate(sites)
        mat .+= _build_local_term(local_terms, local_couplings[i], sitevec, N)
    end

    return mat
end


function _make_diffs(sitevec, N)

    diffs = diff([0, sitevec..., N-1])#[1; sitevec; N]
    diffs[2:end-1] .-= 1
    diffs .= 2 .^diffs 
    return diffs
end

# N = 12
# @benchmark build_local_terms([[sx, sy], [sx, sy]], [[1., 1.] for _ in 1:N], [[i, i%N] for i in 1:N], N)