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

## Running on multiple cores
Use this command in the shell to set the number of processors used in the calculation:
```sh
export  JULIA_NUM_THREADS=${THE_NUMBER_OF_THE_CORES}   
```
It should be set before activating the Julia REPL environment.


## Solver API

### Single electron
#### 1. `eigen_solve(grid_size, potential, cell_length, hbar, mass; sparse)`

Solves the eigenvalue problem for a given Hamiltonian matrix. The space will be devided into sparse grid.

#### Parameters:

- `grid_size::Int`: The size of the grid for the problem.
- `potential::Array`: The potential energy array.
- `cell_length::AbstractFloat`: The length of each cell in the grid.
- `hbar::AbstractFloat`: The reduced Planck's constant.
- `mass::AbstractFloat`: The mass of the particle.
- `sparse::Bool` (optional, default `true`): A flag indicating whether to use sparse matrix representation. (it will save memory, especially when thge grid_size is large)

#### Returns:

- `evals`: The eigenvalues of the Hamiltonian matrix.
- `evecs`: The corresponding eigenvectors of the Hamiltonian matrix.

#### 2. `variation_solve(basis_func, num_orbitals, start_bound, end_bound; potential, symmetric, dimension)`

Find the ground state energy using the variational method.

#### Parameters:

- `basis_func::Function`: The basis function for the problem.
- `num_orbitals::Int`: The number of orbitals in the system.
- `start_bound::Vector`: The starting boundary for the problem.
- `end_bound::Vector`: The ending boundary for the problem.
- `potential::Function` (optional, default `(_) -> 0`): The potential energy function.
- `symmetric::Symbol` (optional, default `:Nothing`): A flag indicating whether the problem is symmetric.
- `dimension::Int` (optional, default `1`): The dimension of the problem.

#### Returns:

- `evals`: The eigenvalues of the Hamiltonian matrix.
- `evecs`: The corresponding eigenvectors of the Hamiltonian matrix.

### Multiple electrons
#### 3. `hartree_fock_solve(nuclear_z, basis_func, num_orbitals; symmetric, dimension, iter_steps, torrlence, start_bound, end_bound, potential_func, mode, converge_alpha)`

Solves the Hartree-Fock equations for a given system.

#### Parameters:

- `nuclear_z::Int`: The atomic number of the nucleus.
- `basis_func::Function`: The basis function for the problem.
- `num_orbitals::Int`: The number of orbitals in the system.
- `symmetric::Symbol` (optional, default `:Spherical`): The symmetry of the system.
- `dimension::Int` (optional, default `3`): The dimension of the problem.
- `iter_steps::Int` (optional, default `50`): The number of iteration steps.
- `torrlence::Float64` (optional, default `1e-2`): The tolerance for convergence.
- `start_bound::Vector` (optional, default `[0.0]`): The starting boundary for the problem.
- `end_bound::Vector` (optional, default `[4.0]`): The ending boundary for the problem.
- `potential_func::Function` (optional, default `(r) -> -nuclear_z / r`): The potential energy function.
- `mode::Symbol` (optional, default `:RHF`): The mode of the calculation (only RHF, that is, restricted Hartree–Fock, is implemented).
- `converge_alpha::Float64` (optional, default `0.6`): The convergence factor asigning the weight of the density matrix in last step and last-last step

#### Returns:

- `ground_state_energy`: The energy of the ground state.


## Example (Test)

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
2. Set the positive charge distribution as a parameter instead of fixing it in the origin
