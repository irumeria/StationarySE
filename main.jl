include("src/StationarySE.jl")

using .StationarySE

println("Runing on ", Threads.nthreads(), " cores")

hydrogen_atom()
