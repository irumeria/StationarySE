using StationarySE
using Test
using LinearAlgebra

@testset "StationarySE.jl" begin

	@info "change cores number by setting JULIA_NUM_THREADS environment variable"
	println("Running on ", Threads.nthreads(), " cores")

	println("")
	println("=== Testing H atom ===")
	potential_func = (r) -> -1/r
	evals, _ = variation_solve(gaussian_hydrogen, 4, [0.0], [3.0]; potential=potential_func, symmetric=:Spherical, dimension=3)
	ground_state_energy = evals[1]
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for H atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")
	println("")

	println("")
	println("=== Testing He atom ===")
	ground_state_energy = hartree_fock_solve(2, gaussian_helium, 4; end_bound=[6.0])
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for He atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")

	println("")
	println("=== Testing C atom ===")
	ground_state_energy = hartree_fock_solve(6, gaussian_carbon, 6, end_bound=[6.0])
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for C atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")

	println("")
	println("=== Testing O atom ===")
	ground_state_energy_1 = hartree_fock_solve(8, gaussian_oxygen_constrain, 8, end_bound=[8.0])
	energy_ev_1 = ground_state_energy_1 * 27.211396
	ground_state_energy_2 = hartree_fock_solve(8, gaussian_oxygen_unconstrain, 8, end_bound=[8.0])
	energy_ev_2 = ground_state_energy_2 * 27.211396
	println("=== different basis functions will give different results:")
	println("(basis 1) O atom: ", ground_state_energy_1, " hatree, or ", energy_ev_1, " ev")
	println("(basis 2) O atom: ", ground_state_energy_2, " hatree, or ", energy_ev_2, " ev")
end
