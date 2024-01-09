using StationarySE
using Test
using GSL
using LinearAlgebra

@testset "StationarySE.jl" begin

	println("testing the Legendre Expansion of potential")
	# normal integral
	f = (s, r, params) -> begin
		dr = norm(r[1:3] - r[4:6])
		s[1] = 1/dr * 1/norm(r[1:3]) * 1/norm(r[4:6])
		s
	end
	res1 = monte_carlo_integral(f, repeat([-10.0], 6), repeat([10.0], 6), 10000000)
	@show res1

	# Legendre Expansion
	f = (s, r, params) -> begin
		lmax = params[1]
		sum = 0.0
		theta = r[3]
		r1 = r[1]
		r2 = r[2]
		rmin = min(r1, r2)
		rmax = max(r1, r2)
		for l in 0:lmax
			# @show theta l
			factor2 = GSL.sf_legendre_Pl(l, cos(theta))
			factor1 = rmin^l / rmax^(l+1)
			sum += factor1 * factor2
		end
		s[1] = sum * 1/r1 * 1/r2
		return s
	end

	res2 = monte_carlo_integral(f, [0, 0, 0], [10, 10, pi], 10000000; params=[100])
	@show res2

    println("Running on ", Threads.nthreads(), " cores")
	nuclear_z = 2
	wave_func = gaussian_helium
	symmetric = :Spherical
	dimension = 3
	iter_steps = 50
	torrlence = 1e-2
	start_bound = [0.0]
	end_bound = [10.0]
	potential_func = (r) -> -2/r
	num_orbitals = 4
	wave_func = gaussian_helium

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

	Q_matrix = two_electron_matrix(
		wave_func,
		num_orbitals,
		start_bound,
		end_bound;
		operator = (dr) -> 1.0/dr,
		dimension = dimension,
		symmetric = symmetric,
		normalize_vec = ec
	)

	@show Q_matrix

	fock_matrix = h_matrix + Q_matrix

	transform_matrix = transform(overlap_matrix)
	evals, evecs = solve_energy(
		change_basis(
		fock_matrix, transform_matrix
		),
	)
	# @show evals
	compute_ground_energy = (evals, Q, coefs, H) -> begin
		2 * coefs' * H * coefs + coefs' * (Q * coefs)
	end

	ground_state_energy = compute_ground_energy(
		evals, Q_matrix, ec, h_matrix)

	@show ground_state_energy

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

		Q_matrix = two_electron_matrix(
			wave_func,
			num_orbitals,
			start_bound,
			end_bound;
			operator = (dr) -> 1.0/dr,
			dimension = dimension,
			symmetric = symmetric,
			normalize_vec = ec
		)

		fock_matrix = h_matrix + Q_matrix
		transform_matrix = transform(overlap_matrix)
		evals, evecs = solve_energy(
			change_basis(
				fock_matrix, 
				transform_matrix
			)
		)
		ground_state_energy = compute_ground_energy(
			evals, Q_matrix, ec, h_matrix)

		println("ground state energy in epoch ", i, " : ", ground_state_energy)

	end
	
	@show ground_state_energy
	
	println("total ground state energy: ", ground_state_energy*2) 
end
