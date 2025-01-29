"""
    interaction(state::T, i::T, j::T) -> Tuple{T, T} where {T <: Integer}

Calculates the interaction term for the given state.

# Arguments
- `state::T`: The quantum state.
- `i::T`: The first site index.
- `j::T`: The second site index.

# Returns
- `Tuple{T, T}`: A tuple containing the product of the factors and the unchanged state.
"""
@inline function interaction(state::T, i::T, j::T) where {T <: Integer}

    factor_i, _ = ops.sz(state, i)
    factor_j, _ = ops.sz(state, j)

    return factor_i * factor_j, state

end


"""
    spsm(state::T, i::T, j::T) -> Tuple{T, T} where {T <: Integer}

Applies the spin raising and lowering operators to the given state.

# Arguments
- `state::T`: The quantum state.
- `i::T`: The site index for the raising operator.
- `j::T`: The site index for the lowering operator.

# Returns
- `Tuple{T, T}`: A tuple containing the product of the factors and the new state.
"""
@inline function spsm(state::T, i::T, j::T) where {T <: Integer}

    factor_sm, state = ops.sm(state, j)
    factor_sp, state = ops.sp(state, i)
    return factor_sm * factor_sp, state
end

"""
    smsp(state::T, i::T, j::T) -> Tuple{T, T} where {T <: Integer}

Applies the spin lowering and raising operators to the given state.

# Arguments
- `state::T`: The quantum state.
- `i::T`: The site index for the lowering operator.
- `j::T`: The site index for the raising operator.

# Returns
- `Tuple{T, T}`: A tuple containing the product of the factors and the new state.
"""
@inline function smsp(state::T, i::T, j::T) where {T <: Integer}

    factor_sp, state = ops.sp(state, j)
    factor_sm, state = ops.sm(state, i)
    return factor_sm * factor_sp, state
end


# Added this for code organization; it appears beneficial to the performance
# if the iterable list is cast as a tuple, not as a vector; this is probably
# due to allocations which can arise if the functions spsm and smsp are
# type-instable; in that case, if the function is called many times inside a loop,
# (which it is in our case), the compiler will have to allocate memory for the return
# values
const _operators = (spsm, smsp)

"""
An internal routine to estimate the number of nonzero matrix elements. We
are quite conservative here.
"""
function _prepare_nnz(ham::XXZHamiltonian{T, N, M}) where {T, N, M}

    nnz = ham.basis.nstates * (2 * ham.basis.size + 1)

    return nnz

end

"""
    _apply_diag!(ham, i, state, site, rows, cols, vals, count) -> Integer

Applies the diagonal terms of the Hamiltonian to the given state.

# Arguments
- `ham`: The XXZHamiltonian instance.
- `i`: The state index.
- `state`: The quantum state.
- `site`: The site index.
- `rows`: The row indices of the Hamiltonian matrix.
- `cols`: The column indices of the Hamiltonian matrix.
- `vals`: The values of the Hamiltonian matrix.
- `count`: The current count of non-zero elements.

# Returns
- `Integer`: The updated count of non-zero elements.
"""
@inline function _apply_diag!(ham, i, state, site, rows, cols, vals, count)

    ifactor, _ = interaction(state, site, (site+1) % ham.basis.size)
    ffactor, _ = ops.sz(state, site)
    count += 1
    rows[count] = i
    cols[count] = i
    vals[count] = ifactor * ham.Delta * ham.J + ffactor * ham.fields[site + 1]

    return count
end

"""
    _apply_offdiag!(ham, op, i, state, site, rows, cols, vals, count) -> Integer

Applies the off-diagonal terms of the Hamiltonian to the given state.

# Arguments
- `ham`: The XXZHamiltonian instance.
- `op`: The operator function to apply.
- `i`: The state index.
- `state`: The quantum state.
- `site`: The site index.
- `rows`: The row indices of the Hamiltonian matrix.
- `cols`: The column indices of the Hamiltonian matrix.
- `vals`: The values of the Hamiltonian matrix.
- `count`: The current count of non-zero elements.

# Returns
- `Integer`: The updated count of non-zero elements.
"""
@inline function _apply_offdiag!(ham, op, i, state, site, rows, cols, vals, count)

    factor, newstate = op(state, site, (site+1) % ham.basis.size)
    # check if the matelt exists
    if !iszero(factor)
        count += 1
        val = factor * ham.J * 0.5
        # a short way of doing if clauses: first is the condition followed by ?
        # the first value if true, the second if false
        # we add +1 to the newstate, again, since julia is 1-based
        _idx = searchsortedfirst(ham.basis.states, newstate)
        # construct lower triangular matrix
        rows[count] = i
        cols[count] = _idx
        vals[count] = val
    end # if clause

    return count

end

"""
    buildham!(ham::XXZHamiltonian{T, N, M}) where {T, N, M}

Builds the Hamiltonian matrix for the given XXZHamiltonian instance.

# Arguments
- `ham::XXZHamiltonian{T, N, M}`: The XXZHamiltonian instance.

# Returns
- `Nothing`: The function modifies the `mat` field of the `ham` instance in-place.

# Example
```julia

using QuantumJulia

basis = Basis(4, 2)
ham = XXZHamiltonian(basis, 1.0, 1.0)
buildham!(ham)
```

# Notes
For better code organization and performance, we break the code down into smaller units. 
Here, we group the Hamiltonian terms into separate functions for diagonal and off-diagonal 
terms. This modular approach helps in maintaining the code and can also improve performance 
by allowing the compiler to better optimize smaller, focused functions.
"""
function buildham!(ham::XXZHamiltonian{T, N, M}) where {T, N, M}

    # this is important! Generally, it's even better if one does not
    # allocate inside functions, but rather pass the arrays as arguments.
    # For the sake of simplicity (and there won't be much of a performance hit)
    # we do it inside this time. However, we still preallocate the space needed
    # for these matrices, since we can estimate the number of the matrix elements.
    # This way, we wont have to push elements to the arrays, which is slow.
    nnz = _prepare_nnz(ham)
    rows = Vector{N}(undef, nnz)
    cols = Vector{N}(undef, nnz)
    vals = Vector{T}(undef, nnz)

    count = 0
    # does not check whether the loop is within bounds or not
    # we use enumerate since states go from 0 onwards and julia
    # is 1-based
    @inbounds for (i::N, state) in enumerate(ham.basis.states)
        @inbounds for site in 0:ham.basis.size-1
            # the interacting part + fields part - we group them together
            count = _apply_diag!(ham, i, state, site, rows, cols, vals, count)

            # the spin flipping part could be nicer, but that'll do for the
            # pedagogic purposes
            for op in _operators
               count = _apply_offdiag!(ham, op, i, state, site, rows, cols, vals, count)
            end # over flipping operators
        end # over sites
    end # over states

    # we do the resizing to the actual matrix size at the end
    resize!(rows, count)
    resize!(cols, count)
    resize!(vals, count)

    # this will create a sparse array; the individual rows, cols, vals won't live outside
    # the scope of this function.
    ham.mat = sparse(rows, cols, vals, ham.basis.nstates, ham.basis.nstates)
    ham._isset = true
end