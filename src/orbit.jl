
function normalize_overlap(overlap_matrix, expansion_coefficients)
	expansion_coefficients /
		(expansion_coefficients' * overlap_matrix * expansion_coefficients)^0.5
end

function _extract_bound(start_bound, end_bound, dimension)
	x_start, x_end, y_start, y_end, z_start, z_end = 0, 0, 0, 0, 0, 0
	if dimension == 1
		x_start = start_bound[1]
		x_end = end_bound[1]
	elseif dimension == 2
		x_start = start_bound[1]
		x_end = end_bound[1]
		y_start = start_bound[2]
		y_end = end_bound[2]
	elseif dimension == 3
		x_start = start_bound[1]
		x_end = end_bound[1]
		y_start = start_bound[2]
		y_end = end_bound[2]
		z_start = start_bound[3]
		z_end = end_bound[3]
	end
	x_start, x_end, y_start, y_end, z_start, z_end
end


function integral_overlap(
	wave_func::Function,
	matrix_size::Int64,
	start_bound::Vector, # x_start, y_start, z_start
	end_bound::Vector; # x_end, y_end, z_end
	operator::Function = (_) -> 1, # Deafult : I
	combined::Bool = false,
	dimension::Int64 = 1,
	symmetric::Symbol = nothing,
	integral_steps::Int64 = 1000000)
	"""
	calculate ∫_x ∫_x <ψ_n|x><x|O|x><x|ψ_m> numerically
	the operator O (Deafult : I) should not contain any differential terms, e.g. d/dx
	For H operator, use `integral_hamiltonian` instead.
	"""
	int_dimension = dimension
	if symmetric == :Spherical
		int_dimension = 1
	end
	
	if combined
		f = wave_func
	elseif symmetric == :Spherical
		f = (y, r, p) -> begin
			y[1] = wave_func(r[1], p[1])' * operator(r[1]) * wave_func(r[1], p[2]) * 4 * pi * r[1]^2
		end
	elseif dimension == 1
		f = (y, r, p) -> begin
			y[1] = wave_func(r[1], p[1])' * operator(r[1]) * wave_func(r[1], p[2])
		end
	elseif dimension == 2
		f = (y, r, p) -> begin
			y[1] = wave_func(r[1], r[2], p[1])' * operator((r[1], r[2])) * wave_func(r[1], r[2], p[2])
		end
	elseif dimension == 3
		f = (y, r, p) -> begin
			y[1] = wave_func(r[1], r[2], r[3], p[1])' * operator((r[1], r[2], r[3])) * wave_func(r[1], r[2], r[3], p[2])
		end
	end

	overlap_matrix = zeros(matrix_size, matrix_size)
	start_bound .+= rand() * 1e-4
	end_bound .+= rand() * 1e-4
	Threads.@threads for n ∈ 0:matrix_size-1
		sol_vec = zeros(matrix_size)
		for m ∈ 0:matrix_size-1

			prob = IntegralProblem(f, start_bound, end_bound, [n, m]) 
			sol = solve(prob, CubatureJLh(), reltol=1e-3, abstol=1e-3)
			sol_vec[m+1] = sol.u[1]
		end
		overlap_matrix[n+1, :] = sol_vec
	end
	overlap_matrix
end

function from_integraled(
	int_wave_func::Function,
	matrix_size::Int64,
)
	"""
	construct overlap matrix from matrix element <ψ_n|ψ_m> or <ψ_n|O|ψ_m>
	"""
	overlap_matrix = zeros(matrix_size, matrix_size)
	for n ∈ 0:matrix_size-1
		for m ∈ 0:matrix_size-1
			overlap_matrix[n+1, m+1] = int_wave_func(n, m)
		end
	end
	overlap_matrix
end

function transform(overlap_matrix::Matrix)
	eval, evec = eigen(overlap_matrix)
	if all(eval .> 0) 
		rev_ev = 1.0 ./ sqrt.(eval)
	else
		@warn "The overlap matrix with element <ψ_m><ψ_n> is not a positive definite matrix."
		rev_ev = 1.0 ./ sqrt.(Complex.(eval))
	end
	diag_matrix = Diagonal(rev_ev)
	evec * diag_matrix
end

export Orbit, gaussian_basis, normalize_overlap, from_gaussian_orbit, integral_overlap, from_integraled, transform
