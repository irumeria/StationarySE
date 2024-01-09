module StationarySE

using LinearAlgebra
using Base.Threads
using Symbolics
using SparseArrays
using KrylovKit
using GeneralizedGenerated
using Integrals
using Cubature
using SharedArrays
using ProgressMeter

include(joinpath(@__DIR__, "math/integral.jl"))

include(joinpath(@__DIR__, "orbit.jl"))
include(joinpath(@__DIR__, "hartreefock.jl"))
include(joinpath(@__DIR__, "hamiltonian.jl"))
include(joinpath(@__DIR__, "solvers.jl"))
include(joinpath(@__DIR__, "buildin.jl"))
include(joinpath(@__DIR__, "demos.jl"))

end
