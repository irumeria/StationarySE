# matrix of <ψ_p|x><ψ_r|x><x|O|x><x|ψ_q><x|ψ_s>
# the |x> integration is done before calling this function
function two_electron_matrix_from_integraled(
	int_fn::Function,
	matrix_size,
	coefs = nothing)

	if isnothing(coefs)
		coefs = ones(matrix_size)
	elseif any(isa(x, Complex) for x in coefs)
		coefs = sqrt.(real.(coefs)^2 + imag.(coefs)^2)
		@warn "coefs contains complex number, only real part is used."
	end

	_mat = zeros((matrix_size, matrix_size, matrix_size, matrix_size))
	mat = zeros((matrix_size, matrix_size))
	for i ∈ 0:matrix_size-1
		for j ∈ 0:matrix_size-1
			for k ∈ 0:matrix_size-1
				for l ∈ 0:matrix_size-1
					_mat[i+1, j+1, k+1, l+1] = int_fn(i, j, k, l)
				end
			end
		end
	end
	
	for i ∈ 0:matrix_size-1
		for j ∈ 0:matrix_size-1
			mat[i+1, j+1] = coefs' * _mat[i+1, :, j+1, :] * coefs
		end
	end
	mat
end

function two_electron_matrix(
	wave_func::Function,
	matrix_size::Int64,
	coefs::Matrix,
	start_bound::Vector, # x_start, y_start, z_start
	end_bound::Vector; # x_end, y_end, z_end
	operator::Function = (_) -> 1, # Deafult : I
	combined::Bool = false,
	dimension::Int64 = 1,
	symmetric::Symbol = nothing)
	"""
	matrix of <ψ_p|x><ψ_m|x><x|O|x><x|ψ_q><x|ψ_s>
	"""
	fock_levels_num = size(coefs)[2]
	for en ∈ 1:fock_levels_num
		if any(isa(x, Complex) for x in coefs[:,en])
			coefs[:,en] = sqrt.(real.(coefs[:,en])^2 + imag.(coefs[:,en])^2)
			@warn "normalize factors contains complex number, the norm is used."
		end
	end

	if combined
		f = wave_func
	elseif symmetric == :Spherical
		f = (s, r, p) -> begin
			r1 = r[1]
			r2 = r[2]
			cos_theta = cos(r[3])
			norm_r2_minus_r1 = sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cos_theta)
			s[1] = wave_func(r1, p[1])' * 
				wave_func(r2, p[2])' * 
				operator(norm_r2_minus_r1) * 
				wave_func(r1, p[3]) * 
				wave_func(r2, p[4]) * 
				4 * pi * r1^2 *
				2 * pi * r2^2 * sin(r[3])
			s
		end
		start_bound = [start_bound[1], start_bound[1], 0]
		end_bound = [end_bound[1], end_bound[1], pi]
	else
		ErrorException("No implement for dimension "+string(dimension))
	end

	shared_mat = SharedArray{Float64}((matrix_size, matrix_size, matrix_size, matrix_size))

	permutations_array = []
	J_mat = zeros((matrix_size, matrix_size))
	K_mat = zeros((matrix_size, matrix_size))
	for i ∈ 0:matrix_size-1
		for j ∈ 0:matrix_size-1
			for k ∈ 0:matrix_size-1
				for l ∈ 0:matrix_size-1
					permutations_array = [permutations_array; [i j k l]]
				end
			end
		end
	end

	@info "Start calculating two electron overlap matrix..."
	prog = Progress(size(permutations_array)[1], barglyphs = BarGlyphs("[=> ]"), barlen = 50, color = :yellow)

	Threads.@threads for p in 1:size(permutations_array)[1]
		i, j, k, l = permutations_array[p, :]
		prob = IntegralProblem(f, start_bound, end_bound, [i, j, k, l]) 
		sol = solve(prob, CubatureJLh(), reltol=1e-3, abstol=1e-3).u[1]
		shared_mat[i+1, j+1, k+1, l+1] = sol
		next!(prog)
	end
	finish!(prog)
	_mat = copy(shared_mat)

	for i ∈ 1:matrix_size
		for j ∈ 1:matrix_size
			for k ∈ 1:fock_levels_num
				J_mat[i, j] += coefs[:, k]' * _mat[i, :, j, :] * coefs[:, k]
			end
		end
	end

	for i ∈ 1:matrix_size
		for j ∈ 1:matrix_size
			for k ∈ 1:fock_levels_num
				K_mat[i, j] += coefs[:, k]' * _mat[i, :, :, j] * coefs[:, k]
			end
		end
	end

	J_mat, K_mat

end

export two_electron_matrix_from_integraled, two_electron_matrix
