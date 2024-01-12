using StationarySE
using Test
using GSL
using LinearAlgebra

@testset "StationarySE.jl" begin
	println("Running on ", Threads.nthreads(), " cores")
	@info "change cores number by setting JULIA_NUM_THREADS environment variable"

	println("")
	println("=== Testing H atom ===")
	potential_func = (r) -> -1/r
	evals, _ = variation_solve(gaussian_hydrogen, 4, [0.0], [100.0]; potential=potential_func, symmetric=:Spherical, dimension=3)
	ground_state_energy = evals[1]
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for H atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")
	println("")

	println("")
	println("=== Testing He atom ===")
	ground_state_energy = hartree_fock_solve(2, gaussian_helium, 2; end_bound=[8.0])
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for He atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")

	println("")
	println("=== Testing C atom ===")
	ground_state_energy = hartree_fock_solve(6, gaussian_carbon, 8; end_bound=[8.0])
	energy_ev = ground_state_energy * 27.211396
	println("total ground state energy for C atom: ", ground_state_energy, " hatree, or ", energy_ev, " ev")
end
