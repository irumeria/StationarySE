
function variation_solve(
	basis_func::Function,
	num_orbitals::Int,
	start_bound::Vector,
	end_bound::Vector;
	potential::Function = (_) -> 0,
	symmetric::Symbol = :Nothing,
	dimension::Int = 1,
)
	overlap_matrix = integral_overlap(
		basis_func,
		num_orbitals,
		start_bound,
		end_bound;
		dimension = dimension,
		symmetric = symmetric,
		)

	h_matrix = integral_hamiltonian(
		basis_func,
		num_orbitals,
		start_bound,
		end_bound;
		potential_func = potential,
		symmetric = symmetric, 
		dimension = dimension 
	)

	new_h = change_basis(h_matrix, transform(overlap_matrix))
	evals, evecs = solve_energy(new_h)
	evals, evecs
end

function eigen_solve(
        grid_size::Int,
        potential::Array,
        cell_length::AbstractFloat,
        hbar::AbstractFloat,
        mass::AbstractFloat;
		sparse = true
    )

	h_matrix = constuct_h_matrix(
    	cell_length, 
        grid_size, 
        hbar, 
        mass, 
        potential;
		sparse = sparse
    )

	evals, evecs = solve_energy(h_matrix)
	evals, evecs
end

function hartree_fock_solve(
		nuclear_z::Int,
		basis_func::Function,
		num_orbitals::Int;
		symmetric = :Spherical,
		dimension = 3,
		iter_steps = 50,
		torrlence = 1e-2,
		start_bound = [0.0],
		end_bound = [4.0],
		potential_func = (r) -> -nuclear_z / r,
		mode = :RHF,
		converge_alpha = 0.6,
	)

	if mode != :RHF
		ErrorException("Only implement RHF mode for now")
	end

	if symmetric != :Spherical
		ErrorException("Only implement spherical symmetry for now, the positive charges are fix in the origin.")
	end

	if nuclear_z % 2 != 0
		ErrorException("RHF mode only supports close-shell system, nuclear_z should be an odd number")
	end

	overlap_matrix = integral_overlap(
		basis_func,
		num_orbitals,
		start_bound,
		end_bound;
		dimension = dimension,
		symmetric = symmetric,
	)
	ec = ones(num_orbitals, nuclear_z ÷ 2)
	for i ∈ 1:size(ec)[2]
		ec[:, i] = normalize_overlap(overlap_matrix, ec[:, i])
	end

	h_matrix = integral_hamiltonian(
		basis_func,
		num_orbitals,
		start_bound,
		end_bound;
		potential_func = potential_func,
		symmetric = symmetric,
		dimension = dimension,
	)
	density_matrix = generate_density_matrix(ec)
	J_matrix, K_matrix = two_electron_matrix(
		basis_func,
		num_orbitals,
		density_matrix,
		start_bound,
		end_bound;
		operator = (dr) -> 1.0 / dr,
		dimension = dimension,
		symmetric = symmetric,
	)

	Q_matrix = 2 * J_matrix - K_matrix

	compute_ground_energy = (Q, coefs, H) -> begin
		eng = 0.0
		for k ∈ 1:size(coefs)[2]
			eng += 2 * coefs[:, k]' * H * coefs[:, k] + coefs[:, k]' * Q * coefs[:, k]
		end
		eng
	end

	ground_state_energy = compute_ground_energy(
		Q_matrix, ec, h_matrix
	)

	@show ground_state_energy

	fock_matrix = h_matrix + Q_matrix

	transform_matrix = transform(overlap_matrix)
	evals, evecs = solve_energy(
		change_basis(
			fock_matrix, transform_matrix,
		),
	)

	last_energy = 0.0
	last_density_matrix = zeros(size(density_matrix))
	for i in 1:iter_steps+1

		deltaE = abs(ground_state_energy - last_energy)
		last_energy = ground_state_energy

		if deltaE < torrlence
			println("converged!")
			break
		elseif i == iter_steps + 1
			println("not converged!")
			break
		end

		evecs = transform_matrix * evecs
		for k ∈ size(ec)[2]
			ec[:, k] = real.(evecs[:, argmin(real.(evals))])
		end

		for k ∈ 1:size(ec)[2]
			ec[:, k] = normalize_overlap(overlap_matrix, ec[:, k])
		end
		_density_matrix = generate_density_matrix(ec)
		density_matrix = converge_alpha .* _density_matrix + (1 - converge_alpha) .* last_density_matrix
		last_density_matrix = density_matrix

		J_matrix, K_matrix = two_electron_matrix(
			basis_func,
			num_orbitals,
			density_matrix,
			start_bound,
			end_bound;
			operator = (dr) -> 1.0 / dr,
			dimension = dimension,
			symmetric = symmetric,
		)

		Q_matrix = 2 * J_matrix - K_matrix
		fock_matrix = h_matrix + Q_matrix
		transform_matrix = transform(overlap_matrix)
		evals, evecs = solve_energy(
			change_basis(
				fock_matrix,
				transform_matrix,
			),
		)
		ground_state_energy = compute_ground_energy(
			Q_matrix, ec, h_matrix,
		)

		println("ground state energy in epoch ", i, " : ", ground_state_energy)

	end
	
	ground_state_energy
	
end

export eigen_solve, variation_solve, hartree_fock_solve
