function harmonic_oscillator(grid_size::Int64; dimension = 2)
	# The center is at the half of the grid
	potential = nothing
	center = grid_size / 2
	if dimension == 1
		potential = zeros(grid_size)
		for i in 1:grid_size
			potential[i] = 0.5 * (i - 1 - center)^2
		end
	elseif dimension == 2
		potential = zeros(grid_size, grid_size)
		for i in 1:grid_size
			for j in 1:grid_size
				potential[i, j] = 0.5 * ((i - 1 - center)^2 + (j - 1 - center)^2)
			end
		end
	elseif dimension == 3
		potential = zeros(grid_size, grid_size, grid_size)
		for i in 1:grid_size
			for j in 1:grid_size
				for k in 1:grid_size
					potential[i, j, k] = 0.5 * ((i - 1 - center)^2 + (j - 1 - center)^2 + (k - 1 - center)^2)
				end
			end
		end
	end
	potential
end

function finite_well(grid_size::Int64; dimension = 1)
	potential = nothing
	center = grid_size / 2
	quarter = grid_size / 4
	if dimension == 1
		potential = zeros(grid_size)
		for i in 1:grid_size
			if abs(i - center) < quarter
				potential[i] = 0.0
			else
				potential[i] = 0.0
			end
		end
	elseif dimension == 2
		potential = zeros(grid_size, grid_size)
		for i in 1:grid_size
			for j in 1:grid_size
				if abs(i - center) < quarter && abs(j - center) < quarter
					potential[i, j] = 0.0
				else
					potential[i, j] = 1.0
				end
			end
		end
	elseif dimension == 3
		potential = zeros(grid_size, grid_size, grid_size)
		for i in 1:grid_size
			for j in 1:grid_size
				for k in 1:grid_size
					if abs(i - center) < quarter && abs(j - center) < quarter && abs(k - center) < quarter
						potential[i, j, k] = 0.0
					else
						potential[i, j, k] = 1.0
					end
				end
			end
		end
	end
	potential
end

function pow_series(x, n)
	x^n * (x - 1) * (x + 1)
end

function analysis_infinit_well(x, n)
	if n % 2 == 0
		cos(n * pi / 2 * x)
	else
		sin(n * pi / 2 * x)
	end
end

function gaussian_hydrogen(r, n)
	"""
	the params of alpha were fitted by 
	R.Ditchfield,W.J.Hehre,andJ.A.Pople, Self-consistentmolecularorbitalmethods.VI.Energy
	optimised Gaussian atomic orbitals, J. Chem. Phys., 52 (1970), 5001-7.

	The result of this is exactly correct, even in a small subspace (num_orbit=4)
	but one can also appoximate the ground state energy of hydrogen by guessing a simple function,
	like the function 'hydrogen_basis_another', and expand it to a large subsapce
	"""
	alpha = 13.007730000000008 - 11.842224870790387*n - 3.435560075773208*n^2 + 5.3563317625429665*n^3 - 1.1241978159793828*n^4
	exp(-alpha * r^2) 
end

function gaussian_helium(r, n)
	alpha = 0.29807299999998527 + 1.646466776632318*n - 0.05681922938143689*n^2 - 1.5927093900343725*n^3 + 0.9475558427835066*n^4
	exp(-alpha * r^2) 
end

function hydrogen_basis_another(r, n)
	"""
	this subsapce will not give a good appoximation unless num_orbitals >= 10
	"""
	alpha = 13/(n+1)
	exp(-alpha * r^2) 
end


export pow_series, gaussian_hydrogen, gaussian_helium, hydrogen_basis_another, harmonic_oscillator, finite_well, hydrogen_atom
