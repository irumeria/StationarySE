
function variation_solve(
	wave_func::Function,
	num_orbitals::Int,
	start_bound::Vector,
	end_bound::Vector;
	potential::Function = (_) -> 0,
	symmetric::Symbol = :Nothing,
	dimension::Int = 1,
)
	overlap_matrix = integral_overlap(
		wave_func,
		num_orbitals,
		start_bound,
		end_bound;
		dimension = dimension,
		symmetric = symmetric,
		)

	h_matrix = integral_hamiltonian(
		wave_func,
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


export eigen_solve, variation_solve
