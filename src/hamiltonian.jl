
function constuct_h_matrix(
	dx::Float64,
	grid_size::Int64,
	hbar::Float64,
	mass::Float64,
	potential::Array{Float64};
	sparse = true)

	dimension = length(size(potential))
	flatten_potential = reshape(potential, length(potential))
	matrix_size = grid_size^dimension
	if sparse
		h_matrix = spzeros(matrix_size, matrix_size)
	else
		h_matrix = zeros(matrix_size, matrix_size)
	end
	quant_minus_2 = -hbar^2 / (2.0 * mass) * (-2.0) / dx^2
	quant_minus_4 = -hbar^2 / (2.0 * mass) * (-4.0) / dx^2
	quant_minus_6 = -hbar^2 / (2.0 * mass) * (-6.0) / dx^2
	grid_size_2 = grid_size^2
	quant_1 = -hbar^2 / (2.0 * mass) * (1.0) / dx^2

	if dimension == 1
		for i in 1:matrix_size
			h_matrix[i, i] = quant_minus_2 + flatten_potential[i]
			if i != 1
				h_matrix[i, i-1] = quant_1
			end
			if i != matrix_size
				h_matrix[i, i+1] = quant_1
			end
		end
	elseif dimension == 2
		for i in 1:matrix_size
			h_matrix[i, i] = quant_minus_4 + flatten_potential[i]
			if i % grid_size != 1
				h_matrix[i, i-1] = quant_1
			end
			if i % grid_size != 0
				h_matrix[i, i+1] = quant_1
			end
			if i > grid_size
				h_matrix[i, i-grid_size] = quant_1
			end
			if i <= matrix_size - grid_size
				h_matrix[i, i+grid_size] = quant_1
			end
		end
	elseif dimension == 3
		for i in 1:matrix_size
			h_matrix[i, i] = quant_minus_6 + flatten_potential[i]
			if i % grid_size_2 != 1 && i % grid_size != 1
				h_matrix[i, i-1] = quant_1
			end
			if i % grid_size_2 != 0 && i % grid_size != 0
				h_matrix[i, i+1] = quant_1
			end
			if i % grid_size_2 > grid_size
				h_matrix[i, i-grid_size] = quant_1
			end
			if i % grid_size_2 <= grid_size_2 - grid_size
				h_matrix[i, i+grid_size] = quant_1
			end
			if i > grid_size_2
				h_matrix[i, i-grid_size_2] = quant_1
			end
			if i <= matrix_size - grid_size_2
				h_matrix[i, i+grid_size_2] = quant_1
			end
		end

	end

	h_matrix
end

function solve_energy(h_matrix; howmany::Int64 = 0)
	# sparse or not
	if howmany <= 0
		if isa(h_matrix, SparseMatrixCSC)
			howmany = size(h_matrix)[1] > 10 ? 10 : size(h_matrix)[1]
		end
	end
	if isa(h_matrix, SparseMatrixCSC)
		evals, evec, _ = eigsolve(h_matrix, howmany, :SR)
	else
		evals, evec = eigen(h_matrix)
	end
	evals, evec
end

function change_basis(h_matrix, transform_matrix)
	transform_matrix' * h_matrix * transform_matrix
end

# <ψ_n|x><x|H|x><x|ψ_m>
function integral_hamiltonian(
	wave_func::Function,
	matrix_size::Int64,
	start_bound::Vector, # x_start, y_start, z_start
	end_bound::Vector; # x_end, y_end, z_end
	potential_func::Function = (_) -> 0, # Deafult : no potential
	symmetric::Symbol = nothing,
	dimension::Int64 = 1)

	v_matrix = integral_overlap(
		wave_func,
		matrix_size,
		start_bound,
		end_bound;
		operator = potential_func,
		symmetric = symmetric,
		dimension = dimension)

	@variables x y z n r

	r = sqrt(x^2 + y^2 + z^2)

	if symmetric == :Spherical
		d2_expr = expand_derivatives(
			(Differential(x)^2)(wave_func(r, n)) +
			(Differential(y)^2)(wave_func(r, n)) +
			(Differential(z)^2)(wave_func(r, n)),
		)
	else
		d2_expr = expand_derivatives(
			sum(Differential(s)^2 for s in [x, y, z][1:dimension])(
				wave_func([x, y, z][1:dimension]..., n),
			),
		)
	end
	d2_phi = mk_function(build_function(d2_expr, [x, y, z][1:dimension]..., n))

	# @show d2_expr
	if symmetric == :Spherical
		combined_wave_func = (y, r, p) -> y[1] = wave_func(r[1], p[1])' * 0.5 * d2_phi(r[1], 0.0, 0.0, p[2]) * 4 * pi * r[1] .^ 2 # why 0.5 ?
	elseif dimension == 1
		combined_wave_func = (y, r, p) -> y[1] = wave_func(r[1], p[1])' * d2_phi(r[1], p[2])
	elseif dimension == 2
		combined_wave_func = (y, r, p) -> y[1] = wave_func(r[1], r[2], p[1])' * d2_phi(r[1], r[2], p[2])
	elseif dimension == 3
		combined_wave_func = (y, r, p) -> y[1] = wave_func(r[1], r[2], r[3], p[1])' * d2_phi(r[1], r[2], r[3], p[2])
	end

	# combined_wave_func = (x, n, m) -> wave_func(x, n) * (wave_func(x + 0.001, m) + wave_func(x - 0.001, m) - 2 * wave_func(x, m)) / 0.001^2

	p_matrix = integral_overlap(
		combined_wave_func,
		matrix_size,
		start_bound,
		end_bound;
		combined = true,
		symmetric = symmetric,
		dimension = dimension)

	-p_matrix + v_matrix

end

export integral_hamiltonian, constuct_h_matrix, solve_energy, change_basis
