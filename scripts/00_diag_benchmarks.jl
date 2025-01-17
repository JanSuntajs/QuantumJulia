# # Diagonalization benchmarks
# This script provides tools for benchmarking
# parallelization capabilities of linear algebra
# libraries used in this project. Specifically,
# we focus on the performance aspects of the
# full exact diagonalization of a dense symmetric
# matrix.
# ## Usage
# 
# To run this script and benchmark the performance of
# different thread numbers, `cd` to the project's root
# folder and issue the following command:
#
# `julia --project=. --threads=auto ./scripts/00_diag_benchmarks.jl <Nmin> <Nmax>`
#
# Where `<Nmin>` and `<Nmax>` specify the minimal and maximal matrix size according to
# the rule `2^N`.
#
# To store the results to disk for a quick visualization, do the following (mind the > sign if you
# wish to output to a file):
#
# `julia --project=. --threads=auto ./scripts/00_diag_benchmarks.jl <Nmin> <Nmax> > <output_file>`
#
# *Note:*  the `--threads=auto` flag is used to enable automatic thread detection. 
using MKL
import LinearAlgebra as la
using BenchmarkTools
using Distributions
using Random
# 
# Below we prepare a random matrix.
#
"""
    prep_rand_mat(N::Int)

A function to prepare a random symmetric
matrix of size 2^N x 2^N for the purpose 
of this benchmark.
"""
function prep_rand_mat(N::Int; dist = Uniform(-1., 1.))

    mat = zeros(2^N, 2^N)

    mat .+= rand(dist, 2^N, 2^N)
    mat .+= transpose(mat)
    mat .*= 0.5

    return mat
    
end
# This is for running the benchmarks.
function run_benchmark(nmin, nmax, max_threads=la.BLAS.get_num_threads())

    threadlist = collect(2 .^ (0:floor(Int, log2(max_threads))))
    if !(max_threads in threadlist)
        push!(threadlist, max_threads)
    end
    println("Starting benchmark! \n")
    for nmat in nmin:nmax

        Random.seed!(1)
        mat = prep_rand_mat(nmat)
        for nthread in threadlist
            println("Matrix size: $(2^nmat) x $(2^nmat). Number of threads: $(nthread) \n")
            la.BLAS.set_num_threads(nthread)
            println(@benchmark la.eigvals($mat))
        end # nthread loop
    end # matrix size loop 
end

function main()

    nmin, nmax = parse.(Int, ARGS)
    #= ----------------------------
    
     PRINTING SYSTEM INFO
    
    --------------------------- =#
    println("Minimal matrix size: $(2^nmin); Maximal matrix size: $(2^nmax). \n")
    println("System info: \n isunix: $(Sys.isunix()) \n islinux: $(Sys.islinux()) \n isapple: $(Sys.isapple()) \n")
    println("CPU info: \n SUMMARY: \n $(Sys.cpu_summary()) \n DETAILED: \n $(Sys.cpu_info()) \n")
    println("MEMORY info: \n TOTAL MEMORY [GB]: $(Int(Sys.total_physical_memory()) / 10^9) \n")
    println("BLAS INFO: \n $(la.BLAS.get_config()) \n")

    run_benchmark(nmin, nmax)


end # main
# 

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end