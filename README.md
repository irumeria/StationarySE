# StationarySE


## Usage
The Library can be tested by running the following command:
```sh
julia           # activate the julia enviroment
]               # activate julia Pkg enviroment
activate .      # activate this project
instantiate     # install and precompile the dependencies for this project
test            # running the build-in test of this project
```

Build-in tests are on the `test\runtest.jl`, they are:
```julia
println("")
println("=== Testing H atom ===")
potential_func = (r) -> -1/r
evals, _ = variation_solve(gaussian_hydrogen, 4, [0.0], [100.0]; potential=potential_func, symmetric=:Spherical, dimension=3)
ground_state_energy = evals[1]


println("")
println("=== Testing He atom ===")
ground_state_energy = hartree_fock_solve(2, gaussian_helium, 2; end_bound=[8.0])

println("")
println("=== Testing C atom ===")
ground_state_energy = hartree_fock_solve(6, gaussian_carbon, 6; end_bound=[8.0])
```

The results of these atom should be:
```sh
-0.4992783898562091 hatree for H atom
-2.852043619183923 hatree for He atom
-38.39474787071248 hatree for C atom
```

## On progress

1. Open-shell system for Hartree-fock Method
2. More universal basis functions
3. Set the positive charge distribution as a parameter instead of fixing it in the origin
