
function harmonic_oscillator_2d()

	grid_size = 10
	potential = harmonic_oscillator(grid_size)
	x = range(-1, 1, grid_size) |> collect
	dx = x[2] - x[1]
	hbar = 1e-3
	mass = 1e-4

	evals, _ = eigen_solve(
		grid_size,
		potential,
		dx,
		hbar,
		mass,
	)

	println("energies: ", evals)

end

function variation_finite_well()
	num_orbitals = 10
	int_wave_func(i, j) = begin
		if i + j == 0 || (i + j) % 2 == 0
			2.0 / (i + j + 5.0) -
			4.0 / (i + j + 3.0) +
			2.0 / (i + j + 1.0)
		else
			0.0
		end
	end
	int_hamiltonian(i, j) = begin
		if (i + j) == 0 || (i + j) % 2 == 0
			-8.0 * (1.0 - i - j - 2.0 * i * j) /
			((i + j + 3.0) * (i + j + 1.0) * (i + j - 1.0))
		else
			0.0
		end
	end
	overlap_matrix = from_integraled(
		int_wave_func,
		num_orbitals)
	# @show overlap_matrix
	h_matrix = from_integraled(int_hamiltonian, num_orbitals)

	new_h = change_basis(h_matrix, transform(overlap_matrix))
	_evals, _ = solve_energy(new_h)
	println("analysis energies: ", _evals)

	wave_func = pow_series

	evals, _ = variation_solve(wave_func, num_orbitals, [-1.0], [1.0]; potential=(_)->0, symmetric=:Nothing, dimension=1)
	println("energies: ", evals)

	# @assert all(evals .≈ _evals)
end

function hydrogen_atom()
	num_orbitals = 4
	wave_func = gaussian_hydrogen
	potential_func = (r) -> -1/r
	evals, _ = variation_solve(wave_func, num_orbitals, [0.0], [100.0]; potential=potential_func, symmetric=:Spherical, dimension=3)
	println("energies: ", evals)
end

function helium_ground_state()

	num_orbitals = 4
	alphas = [0.298073, 1.242567, 5.782948, 38.474970] # TODO
	nuclear_z = 2
	wave_func = gaussian_helium
	symmetric = :Spherical
	dimension = 3
	iter_steps = 50
	torrlence = 1e-2
	start_bound = [0.0]
	end_bound = [50.0]
	potential_func = (r) -> -2/r

	# ==========

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

	@show overlap_matrix

	int_T = begin
		(i, j) ->
			3 * alphas[i+1] * alphas[j+1] * pi^1.5 /
			(alphas[i+1] + alphas[j+1])^2.5
	end
	int_V = begin
		(i, j) ->
			-2 * nuclear_z * pi / (alphas[i+1] + alphas[j+1])
	end
	int_H = begin
		(i, j) ->
			int_T(i, j) + int_V(i, j)
	end

	h_matrix = from_integraled(int_H, num_orbitals)

	@show h_matrix

	# Construct Q Matrix 
	int_q_fn = begin
		(i, j, k, l) -> (2 * pi^(2.5) /
						 ((alphas[i+1] + alphas[k+1]) *
						  (alphas[j+1] + alphas[l+1]) *
						  sqrt(alphas[i+1] + alphas[j+1] + alphas[k+1] + alphas[l+1]))
		)
	end
	Q_matrix = two_electron_matrix_from_integraled(
		int_q_fn,
		num_orbitals,
		ec)

	@show Q_matrix

	pesudo_matrix = h_matrix + Q_matrix
	pesudo_matrix

	transform_matrix = transform(overlap_matrix)
	evals, evecs = solve_energy(
		change_basis(
		pesudo_matrix, transform_matrix
		),
	)
	# @show evals
	compute_ground_energy = (evals, Q, coefs, H) -> begin
		2 * coefs' * H * coefs + coefs' * (Q * coefs)
	end

	ground_state_energy = compute_ground_energy(
		evals, Q_matrix, ec, h_matrix)
	# @show ground_state_energy
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
		ec = evecs[:, argmin(real(evals))] |> copy
		
		# @show orbits.expansion_coefficients
		
		ec = normalize_overlap(overlap_matrix, ec)
		@show ec
		Q_matrix = two_electron_matrix_from_integraled(
			int_q_fn,
			num_orbitals,
			ec)

		pesudo_matrix = h_matrix + Q_matrix
		transform_matrix = transform(overlap_matrix)
		evals, evecs = solve_energy(
			change_basis(
				pesudo_matrix, 
				transform_matrix
			)
		)
		ground_state_energy = compute_ground_energy(
			evals, Q_matrix, ec, h_matrix)

	end
	
	@show ground_state_energy

	# =========
	num_orbitals = 8
	wave_func = (r, n) -> exp(-13.0/(n+1) * r^2)

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
		ec,
		start_bound,
		end_bound;
		operator = (dr) -> 1.0/dr,
		dimension = dimension,
		symmetric = symmetric
	)

	@show Q_matrix

	pesudo_matrix = h_matrix + Q_matrix

	transform_matrix = transform(overlap_matrix)
	evals, evecs = solve_energy(
		change_basis(
		pesudo_matrix, transform_matrix
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
			ec,
			start_bound,
			end_bound;
			operator = (dr) -> 1.0/dr,
			dimension = dimension,
			symmetric = symmetric,
		)

		pesudo_matrix = h_matrix + Q_matrix
		transform_matrix = transform(overlap_matrix)
		evals, evecs = solve_energy(
			change_basis(
				pesudo_matrix, 
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

export harmonic_oscillator_2d, variation_finite_well, hydrogen_atom, helium_ground_state
