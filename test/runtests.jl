using StationarySE
using Test
using GSL
using LinearAlgebra

@testset "StationarySE.jl" begin

    println("Running on ", Threads.nthreads(), " cores")
	nuclear_z = 4
	symmetric = :Spherical
	dimension = 3
	iter_steps = 50
	torrlence = 1e-2
	start_bound = [0.0]
	end_bound = [10.0]
	potential_func = (r) -> -nuclear_z/r
	num_orbitals = 4
	wave_func = gaussian_helium
	# num_orbitals = 6
	# wave_func = (r, n) -> - 2 * cos((n+1)*r)

	if nuclear_z%2 != 0
		ErrorException("Only implement for RHF now, nuclear_z should be odd number")
	end

	overlap_matrix = integral_overlap(
		wave_func,
		num_orbitals,
		start_bound,
		end_bound;
		dimension = dimension,
		symmetric = symmetric,
	)
	ec = ones(num_orbitals)
	ec = normalize_overlap(overlap_matrix, ec)

	h_matrix = integral_hamiltonian(
		wave_func,
		num_orbitals,
		start_bound,
		end_bound;
		potential_func = potential_func,
		symmetric = symmetric, 
		dimension = dimension 
	)

	@show h_matrix

	J_matrix, K_matrix = two_electron_matrix(
		wave_func,
		num_orbitals,
		start_bound,
		end_bound;
		operator = (dr) -> 1.0/dr,
		dimension = dimension,
		symmetric = symmetric,
		normalize_vec = ec
	)

	Q_matrix = 2*J_matrix - K_matrix
	@show J_matrix, K_matrix
	@show Q_matrix

	# @show evals
	compute_ground_energy = (Q, coefs, H, nuclear_z) -> begin
		nuclear_z*(coefs' * H * coefs + 1/2 * coefs' * Q * coefs)
	end

	ground_state_energy = compute_ground_energy(
		Q_matrix, ec, h_matrix, nuclear_z)

	@show ground_state_energy

	fock_matrix = h_matrix + Q_matrix .* (nuclear_z/2)

	transform_matrix = transform(overlap_matrix)
	evals, evecs = solve_energy(
		change_basis(
		fock_matrix, transform_matrix
		),
	)

	last_energy = 0.0

	for i in 1:iter_steps+1

		deltaE = abs(ground_state_energy - last_energy)
		last_energy = ground_state_energy

		if deltaE < torrlence
			println("converged!")
			break
		elseif i == iter_steps+1
			println("not converged!")
			break
		end
		
		evecs = transform_matrix * evecs
		ec = evecs[:, argmin(real.(evals))] |> copy
		
		normalize_overlap(overlap_matrix, ec)

		J_matrix, K_matrix = two_electron_matrix(
			wave_func,
			num_orbitals,
			start_bound,
			end_bound;
			operator = (dr) -> 1.0/dr,
			dimension = dimension,
			symmetric = symmetric,
			normalize_vec = ec
		)

		Q_matrix = 2*J_matrix - K_matrix
		fock_matrix = h_matrix + Q_matrix .* (nuclear_z/2)
		transform_matrix = transform(overlap_matrix)
		evals, evecs = solve_energy(
			change_basis(
				fock_matrix, 
				transform_matrix
			)
		)
		ground_state_energy = compute_ground_energy(
			Q_matrix, ec, h_matrix, nuclear_z)

		println("ground state energy in epoch ", i, " : ", ground_state_energy)

	end
	@show nuclear_z*(ec' * h_matrix * ec)
	@show ec' * Q_matrix * ec
	@show ground_state_energy
	@show evals
	# @show nuclear_z*(evals[1]-1/2*ec'*Q_matrix*ec)

	energy_ev = ground_state_energy * 27.211396
		
	println("total ground state energy: ", ground_state_energy, " hatree, or ", energy_ev, " ev") 
end
