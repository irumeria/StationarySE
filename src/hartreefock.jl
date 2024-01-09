# matrix of <ψ_p|x><ψ_r|x><x|O|x><x|ψ_q><x|ψ_s>
# the |x> integration is done before calling this function
function two_electron_matrix_from_integraled(
	int_fn::Function,
	matrix_size,
	normalize_vec = nothing)

	if isnothing(normalize_vec)
		normalize_vec = ones(matrix_size)
	elseif any(isa(x, Complex) for x in normalize_vec)
		normalize_vec = real.(normalize_vec)
		@warn "normalize_vec contains complex number, only real part is used."
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
			mat[i+1, j+1] = normalize_vec' * _mat[i+1, :, j+1, :] * normalize_vec
		end
	end
	mat
end


function two_electron_matrix(
	wave_func::Function,
	matrix_size::Int64,
	start_bound::Vector, # x_start, y_start, z_start
	end_bound::Vector; # x_end, y_end, z_end
	operator::Function = (_) -> 1, # Deafult : I
	combined::Bool = false,
	dimension::Int64 = 1,
	symmetric::Symbol = nothing,
	integral_steps::Int64 = 1000000,
	normalize_vec = nothing)
	"""
	matrix of <ψ_p|x><ψ_m|x><x|O|x><x|ψ_q><x|ψ_s>
	"""
	if isnothing(normalize_vec)
		normalize_vec = ones(matrix_size)
	elseif any(isa(x, Complex) for x in normalize_vec)
		normalize_vec = real.(normalize_vec)
		@warn "normalize factors contains complex number, only the real part is used."
	end

	if combined
		f = wave_func
	elseif symmetric == :Spherical
		f = (s, r, p) -> begin
			r = reshape(r, (3, 2))
			r1 = norm(r[:, 1])
			r2 = norm(r[:, 2])
			dr = norm(r[:, 1] - r[:, 2])
			s[1] = wave_func(r1, p[1])' * wave_func(r2, p[2])' * operator(dr) * wave_func(r1, p[3]) * wave_func(r2, p[4])
			s
		end
		start_bound = repeat(-sqrt.(end_bound), 6) # for test
		end_bound = repeat(sqrt.(end_bound), 6)
	else
		ErrorException("Not implemented for dimension "+string(dimension))
	end

	# start_bound .+= rand() * 1e-5
	# end_bound .+= rand() * 1e-5

	shared_mat = SharedArray{Float64}((matrix_size, matrix_size, matrix_size, matrix_size))

	permutations_array = []
	mat = zeros((matrix_size, matrix_size))
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
		sol = monte_carlo_integral(
			f, 
			start_bound |> copy, 
			end_bound |> copy, 
			integral_steps |> copy; 
			params = [i, j, k, l], 
			save_memory = true)
		shared_mat[i+1, j+1, k+1, l+1] = sol
		next!(prog)
	end
	finish!(prog)
	_mat = copy(shared_mat)

	for i ∈ 0:matrix_size-1
		for j ∈ 0:matrix_size-1
			mat[i+1, j+1] = normalize_vec' * _mat[i+1, :, j+1, :] * normalize_vec
		end
	end

	mat

end

export two_electron_matrix_from_integraled, two_electron_matrix
