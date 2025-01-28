# # Linear algebra operations in Julia
#
# This is a brief overview of the linear algebra functionalities that are readily available to Julia inside the
# [`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) library. In Julia (as in much of scientific computation), dense linear algebra operations are based on [LAPACK](https://www.netlib.org/lapack/), which is itself built on top of basic linear algebra routines known as [BLAS](https://www.netlib.org/blas/). By default, `LinearAlgebra` ships with `OpenBLAS`, but other, more performant options are also available. See for instance [MKL for Julia](https://github.com/JuliaLinearAlgebra/MKL.jl) (which may or may not perform on Apple silicon processors).
#
# Interfaces for working with both `BLAS` and `LAPACK` routines directly are readily available:
#
# * [BLAS functions](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#BLAS-functions)
# * [LAPACK functions](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg-lapack-functions)

# For working with sparse arrays, one typically uses the [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) library. We also note a [host of external packages](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Noteworthy-External-Sparse-Packages), such as the [KrylovKit](https://github.com/Jutho/KrylovKit.jl)
# which we also employ in this project.
#
# Below, we will cover:
#
# * some basic notes on array creation
# * the importance of preallocation
# * vectorization by means of dotted operations
# * the use of views vs. slices
# * (implicitly, throughout the examples) the importance of benchmarking
#
# For our use cases, the above bulletpoints have proven the most important when writing efficient Julia code.

using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Random

# ## Some basic notes

# First things first, let us briefly discuss basic data containers we would 
# typically use in our numerical calculations - vectors, arrays, matrices 
# (which are actually just an alias for 2D arrays), vectors of vectors ... 
# [This](https://www.matecdev.com/posts/julia-array-initialize.html) is a 
# good starting resource for those wishing to know more.

#= to construct a vector: note the square brackets and commas =#
v1 = [1, 2, 3]
println("Type of v1: ", typeof(v1))
v2 = [1., 2, 3]
println("Type of v2: ", typeof(v2))
#= this would also create a vector (of vectors) =#
v3 = [[1, 2, 3], ["s", "d"], [1.]]
println("Type of v3: ", typeof(v3))

# To create an array/ matrix, do the following:

A = [1 2 3; 4 5 6; 7 8 9]
#

#= matrix is a subtype of the more general Array type =#
println(Matrix <: Array)
#= What about vectors? =#
println(Vector <: Array)
#= what about vectors and matrices? =#
println(Vector <: Matrix)

# During construction, bear this in mind: 
#
# when constructing vectors, each element can be separated either 
# by commas or semicolons; however, when separating by whitespace, 
# that creates a matrix, which is a different type of entity in Julia.
#

b1 = [1, 2, 3]
b2 = [1; 2; 3]
b3 = [1 2 3]
println("Type of b1: ", typeof(b1))
println("Type of b2: ", typeof(b2))
println("Type of b3: ", typeof(b3))

# ## Array comprehensions, undefined arrays and array preallocation
#
# Of course there are many different ways to construct arrays other than 
# explicit construction, for instance array comprehensions:

#= the semicolon at the end suppresses the output
we are creating a vector of inverse squares of integers from 1 to 1000 =#
vec = [1/i^2 for i in 1:1000];
#= let's sum this up =#
sum(vec)

# If parentheses are used instead of square brackets, a *generator* object 
# is created. This may sometimes help with the performance, as a generator 
# is not evaluated untill needed. 

genvec = (1/i^2 for i in 1:1000)
sum(genvec)

# Let's benchmark the two approaches:
@btime sum([1/i^2 for i in 1:10000])
#
@btime sum((1/i^2 for i in 1:10000))

# While the former method is faster, it also allocates more. In general, 
# it is advisable to benchmark you code snippets this way to see what yields optimal performance.
# 
# For performance reasons, it is typically recommended to initialize arrays of a given type (and size) 
# without any values (not even zeros). To do this, we use keywords such 
# as `Vector{T}, Matrix{T}, Array{T}`, where `T` is a given type, such as `Float64` or `Int64`. 
# The syntax is as follows:

n = 10
A0 = Array{Float64}(undef, n, n);
A1 = Array{Float64, 2}(undef, n, n);
A2 = Matrix{Int64}(undef, n, n);

V0 = Vector{Float64}(undef, n);

# The above approach will reserve a certain portion of memory for the containers in advance. 
# As it is not initializing them, this saves some time. But beware, as it may also lead to unwelcome 
# surprises later down the road if you assume the existence of some sensible values (in short: do not 
# assume they were initialized to zero, they may just as well be some gibberish, this can 
# be system and compiler dependent).
#
# To show where this is useful, consider filling a vector with a known number of elements, 
# first using preallocation, then without it. 

n = 100000
@btime begin
   #= apparently, using eachindex is also important
   as it makes a difference compared to 1:n looping =#
   @inbounds for i in eachindex(V1)
         V1[i] = i
   end
end setup=(V1 = Vector{Int64}(undef, n))
#
#= now without the preallocation =#
V2 = Int64[]
@btime begin
    for i in 1:n
        push!(V2, i)
    end
end setup=(V2 = Int64[])

# That appears as a huge difference, both in time and memory. Again, it is 
# important to benchmark your code snippets, as the outcomes are sometimes surprising 
# at a first glance (to be honest, sometimes even after a prolonged stare). 
# In the above example, the slowdown was expected, as the vector needs to be 
# reallocated upon each addition of an element. In terms of performance, 
# try to never allocate new collections inside for loops; instead, 
# do preallocation in advance. If possible, also try to reduce the number 
# of allocations inside functions, especially if those are used inside long loops. 

# Finally, there is this little trick which we will keep coming back to:

#= case 1 =#
n = 100000
#= note the dot syntax =#
@btime V3 .= 1:n setup=(V3 = Vector{Int64}(undef, n));

# We will return back to Julia's 
# [dot syntax for vectorizing functions](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized) 
# shortly. Note that to initialize arrays to zero values, one can use:

#= as said before, sometimes the compiler will
also initialize undef to zeros, but it's best
not to expect it to avoid headaches. =#
zeroarr = zeros(Float64, 10, 10)

# ## Vectorization via the dot syntax
#
# *Any* Julia function can be applied elementwise to any array (or other collection)
# using the syntax `f.(A)`. Following the official docs:

A = [1. 2. 3.; 4. 5. 6.; 7. 8. 9.]

sin.(A)

# Consult the docs for a more thorough explanation, we'll make do with some examples here:
f(x, y) = 3x + 4y;

A = [1., 2., 3.]
B = [4., 5., 6.]

r1 = f.(A, B);
r2 = f.(π, B);

# If possible, multiple dotted operations will be *fused* into a single one. In plain terms, 
# this is beneficial, as allocations of intermediate arrays are usually skipped. One can also 
# imagine that the for loops involved in the calculation (that are done under the hood) 
# get merged into a single loop as our final benchmarking example below suggests.
y(x) = x^2
h(x) = cos(x)

y.(h.(A))
#= this would achieve the same as =#
broadcast(x -> y(h(x)), A) == y.(h.(A)) == [y(h(a)) for a in A]

# Benchmarks below are to prove these methods are comparable/the same in terms of performance. 
# This also goes on to show that quite often, vectorization isn't even needed in Julia, as "for loops are fast" 
# (with some caveats).

@btime y.(h.(X)) setup=(X=collect(1:10000));
#
@btime broadcast(x -> y(h(x)), X) setup = (X=collect(1:10000));
#
@btime [y(h(x)) for x in X] setup=(X=collect(1:10000));

# ## Some more notes on fusing (as per official documentation)

# Compare the functions below. Both achieve the same goal, however, one is considerably 
# faster when applied to an array. Which one? Why?


f(x) = 3x.^2 + 4x + 7x.^3;

#= equivalent to 3 .* x.^2 .+ 4 .* x .+ 7 .* x.^3 =#
fdot(x) = @. 3x^2 + 4x + 7x^3; 
#
@btime f(X) setup = (X=collect(1:10000));
#
@btime fdot(X) setup=(X=collect(1:10000));
#
@btime f.(X) setup = (X=collect(1:10000));

# ## Views, slices, to copy or not to copy?

# In Julia, slicing arrays, such as `arr[1:5, :]`, creates a copy of the data being 
# sliced (except on the left-hand side of an assignment, where `arr[1:5, :] = ...`) 
# assigns in-place to that part of an array. 
# [As per official documentation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views), 
# this may be beneficial for performance in case where many operations are 
# performed on the slice (as the pros of working with a smaller contiguous copy rather 
# than with the original large array outweigh the cons of copying). However, this 
# may not be the case if just some simple operations are performed on the slice and 
# hence the costs of allocation and copying can be considerable. 

# An alternative is to use a *view* of the array (in fact, a `SubArray` object) that 
# references the portion of the original array in-place. Writing to a view also 
# modifies the original array in-place. This is achieved by using the `@view` macro 
# on an individual slice or for the whole expression by using `@views`. An example from the 
# official documentation:

#= copies the data from x by creating a slice =#
fcopy(x) = sum(x[2:end-1]);

#= all slices in the expression are considered as views =#
@views fview1(x) = sum(x[2:end-1]);

#= single view, placed directly before the slice =#
fview2(x) = sum(@view x[2:end-1]);

#= random vector with 1e6 elements =#
x = rand(10^6);
#
@btime fcopy(x);
#
@btime fview1(x);
#
@btime fview2(x);
#
# Note the approximately three times faster execution and considerably smaller memory usage 
# when using views. In most cases of interest to us, it is therefore beneficial to use views, 
# but do ensure it actually brings expected performance gains by properly benchmarking!
#
# ## One final note: 

# Use in-place updates rather than out-of-place, if possible.

#= In this function, we add to Z in-place;
note there is no return specified and
how Z will be changed outside the function =#
function inplace_add!(Z, Y; n=100)
    @inbounds for i in 1:n
            Z .+= Y
    end
end
#
#= in this function, on every iteration, a new
array Z is created. At the end, it is returned.
This function does not change its arguments. =#
function outplace_add(Z, Y; n=100)
    @inbounds for i in 1:n
            Z = Z .+ Y
    end

    return Z
end
#
X = rand(1000, 1000)
Y = rand(1000, 1000)
Z1 = copy(X);
Z2 = copy(X);
Z1
#
inplace_add!(Z1, Y)
Z1 #= Z1 is changed =#
#
outplace_add(Z2, Y) #= Z2 is not changed =#
Z2
#
Z2 = copy(X)
Z2 = outplace_add(Z2, Y) #= Z2 is changed =#

Z1 - Z2
#
#= benchmarks of these routines =#
@btime inplace_add!(Z1, Y) setup = (Z1 = copy(X), Y = Y);
#
@btime Z2 = outplace_add(Z2, Y) setup = (Z2 = copy(X));