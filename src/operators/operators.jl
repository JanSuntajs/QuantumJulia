module Operators

    using Distributions: binomial
    using LinearAlgebra

    include("./bitmanip.jl")
    include("./spinops.jl")


    export getbit, bitflip, countsetbits, selectstates
    export sx, sy, sz, id2, sp, sm

end