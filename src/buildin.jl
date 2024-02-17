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
	the params of alpha were suggested by 
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
	"""
	the params of alpha were suggested by 
	R.Ditchfield,W.J.Hehre,andJ.A.Pople, Self-consistentmolecularorbitalmethods.VI.Energy
	optimised Gaussian atomic orbitals, J. Chem. Phys., 52 (1970), 5001-7.
	"""
	alpha = 0.29807299999998527 + 1.646466776632318*n - 0.05681922938143689*n^2 - 1.5927093900343725*n^3 + 0.9475558427835066*n^4
	exp(-alpha * r^2) 
end

function hydrogen_basis_another(r, n)
	"""
	this subsapce will not give a good appoximation unless num_orbitals ~ 9
	"""
	alpha = 13/(n+1)
	exp(-alpha * r^2) 
end


function gaussian_carbon(r, n)
	"""
	the params of alpha were suggested by 
	R.Ditchfield,W.J.Hehre,andJ.A.Pople, Self-consistentmolecularorbitalmethods.VI.Energy
	optimised Gaussian atomic orbitals, J. Chem. Phys., 52 (1970), 5001-7.
	"""
	x = n
	alpha = 446.62209999992734 - 598.006694158179*x + 119.8768753734617*x^2 + 238.060660375486*x^3 - 194.99390715856327*x^4 + 66.59740985677443*x^5 - 11.854795762057295*x^6 + 1.0781859016858948*x^7 - 0.0395744285359738*x^8
	exp(-alpha * r^2) 
end

function gaussian_oxygen_unconstrain(r, n)
	x = n
	alpha = 282.4431 - 605.9857971179763*x + 655.5668146075697*x^2 - 418.2070589896529*x^3 + 161.77781207998268*x^4 - 38.01963633881945*x^5 + 5.277294957534723*x^6 - 0.3961440535515874*x^7 + 0.012354854913194447*x^8
	exp(-alpha * r^2)
end

function gaussian_oxygen_constrain(r, n)
	x = n
	alpha = 816.2436 - 1484.511751554762*x + 1233.0987376660914*x^2 - 582.5450861416666*x^3 + 167.60307820440968*x^4 - 29.878573704166662*x^5 + 3.222681174930555*x^6 - 0.19268834940476187*x^7 + 0.00490270456845238*x^8
	exp(-alpha * r^2)
end

function build_basis_benzene()
	# TODO: the basis of benzene is not correct

	# (1.40272, 0.0, 0.0),
	# (0.70136, 1.21479, 0.0),
	# (-0.70136, 1.21479, 0.0),
	# (-1.40272, 0.0, 0.0),
	# (-0.70136, -1.21479, 0.0),
	# (0.70136, -1.21479, 0.0),
	@variables rv
	cws = (rv) -> begin # carbon basis set
		ec_carbon = [-0.018694, 0.068236, 0.339788, 0.919272, -0.174387, 0.063763]
		location_C_atoms = [1.40272 0.0 0.0; 
							0.70136 1.21479 0.0; 
							-0.70136 1.21479 0.0; 
							-1.40272 0.0 0.0; 
							-0.70136 -1.21479 0.0; 
							0.70136 -1.21479 0.0]
		carbon_alpha = [446.6221, 67.34026, 14.98591, 3.9944722, 1.967742e-2, 1.35529e-1]
		[sum([ec_carbon[n]*exp(-carbon_alpha[n] * norm(rv-location_C_atoms[m,:])^2) for n in 1:6]) for m in 1:6]
	end

	@variables x y z
	print(expand_derivatives(cws([x, y, z])[1]))
	# construct the basis of benzene based on group theory
	# take benzene as a C6 group for simplicity
	benzene_basis_f_ = (rv)-> begin
		benzene_basis = []
		omega = 2 * pi * im / 6
		push!(benzene_basis, (rv)->cws(rv)[1] + cws(rv)[2] + cws(rv)[3] + cws(rv)[4] + cws(rv)[5] + cws(rv)[6])
		push!(benzene_basis, (rv)->cws(rv)[1] - cws(rv)[2] + cws(rv)[3] - cws(rv)[4] + cws(rv)[5] - cws(rv)[6])
		push!(benzene_basis, (rv)->cws(rv)[1] + omega * cws(rv)[2] + omega^2 * cws(rv)[3] + omega^3 * cws(rv)[4] + 
					omega^4 * cws(rv)[5] + omega^5 * cws(rv)[6])
		push!(benzene_basis, (rv)->cws(rv)[1] + omega^5 * cws(rv)[2] + omega^4 * cws(rv)[3] + omega^3 * cws(rv)[4] + 
					omega^2 * cws(rv)[5] + omega * cws(rv)[6])
		push!(benzene_basis, (rv)->cws(rv)[1] + omega^2 * cws(rv)[2] + omega^4 * cws(rv)[3] + cws(rv)[4] + 
					omega^2 * cws(rv)[5] + omega^4 * cws(rv)[6])
		push!(benzene_basis, (rv)->cws(rv)[1] + omega^4 * cws(rv)[2] + omega^2 * cws(rv)[3] + cws(rv)[4] + 
					omega^4 * cws(rv)[5] + omega^2 * cws(rv)[6])
		[benzene_basis[i](rv) for i in 1:6]
	end
	
	@variables x y z n
	benzene_basis_f_expr = expand_derivatives(benzene_basis_f_([x, y, z])[n])
	benzene_basis_f = mk_function(build_function(benzene_basis_f_expr, [x, y, z, n]))
	print(benzene_basis_f_expr)
	benzene_basis_f

end

export pow_series, gaussian_hydrogen, gaussian_helium, hydrogen_basis_another, harmonic_oscillator, finite_well, hydrogen_atom, gaussian_carbon, gaussian_oxygen_unconstrain, gaussian_oxygen_constrain, build_basis_benzene
