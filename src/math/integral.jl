
function simpson(f, a, b, n; params = [])
	h = (b - a) / n
	x = range(a, stop = b, length = n + 1)
	y = f.(x, params...)
	# when meet NaN in i, use y[i] = (y[i+1] + y[i-1])/2 instead
	# TODO: use a better way to deal with NaN
	for i in 1:n
		if isnan(y[i])
			if i == 1
				y[i] = y[i+1]
			elseif i == n
				y[i] = y[i-1]
			else
				y[i] = (y[i+1] + y[i-1]) / 2
			end
		end
	end
	return h / 3 * (y[1] + 4 * sum(y[2:2:end-1]) + 2 * sum(y[3:2:end-2]) + y[end])
end

function simpson2(f, be, ed, m, n; params = [])
	a, b = be
	c, d = ed
	hx = (b - a) / m
	hy = (d - c) / n
	x = range(a, stop = b, length = m + 1)
	y = range(c, stop = d, length = n + 1)
	integral_sum = 0.0
	
	for i in 1:m
		for j in 1:n
			xi = x[i]
			xj = y[j]
			fij = f(xi, xj, params...)
			
			if isnan(fij)
				if i == 1
					fij = f(x[i+1], xj, params...)
				elseif i == m
					fij = f(x[i-1], xj, params...)
				else
					fij = (f(x[i+1], xj, params...) + f(x[i-1], xj, params...)) / 2
				end
			end
			
			integral_sum += fij * ((hx * hy) / 9) * (4 - (-1)^(i+j))
		end
	end
	
	return integral_sum
end

function monte_carlo_integral(f, be, ed, n; params = [], save_memory = false, in_sphere = false)
	"""
	get dimension from be
	"""
	if !save_memory
		sample_points = rand(Float64, (size(be)[1]), n) .* (ed - be) .+ be
	end
	integral_sum = 0.0
	final_n = n
	
	for i in 1:n

		if save_memory
			sp = rand(Float64, (size(be)[1])) .* (ed - be) + be
		else
			sp = sample_points[:, i]
		end

		if in_sphere
			if norm(sp) > ed[1] - be[1]
				final_n -= 1
				continue
			end
		end
		
		fi = f([0.0], sp, params)[1]
		
		if isnan(fi)
			final_n -= 1
			continue
		end
		
		integral_sum += fi
	end
	volumn = prod(ed - be)
	(integral_sum / final_n) * volumn
end

export simpson, simpson2, monte_carlo_integral
