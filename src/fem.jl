function diffusion_fem(
	s::OVector64,
	pj::Int64,
	ω::Function,
	R::Function,
	D::Function,
	u∞::Float64)

	# Spatial Grid
	P = length(s) - 1
	sj = s[pj]
	# ℙ = 0:P

	# function spatial_grid(p::Integer)
	# 	p ∉ ℙ && error("Index outside grid.")
	# 	p == 0 && return 0.0
	# 	p == pj && return sj
	# 	p == P && return sP
	# 	p < pj && return sj * (1 - (1 - p / pj)^3)
	# 	p > pj && return sj + (sP - sj) * ((p - pj) / (P - pj))^3
	# 	return NaN
	# end

	# s = OffsetArray(
	# 	[spatial_grid(p) for p ∈ ℙ],
	# 	Origin(0)
	# )
	h = s |> parent |> diff

	# Mass & Stiffness Matrices
	R′(s) = R(s) * sin(ω(s))
	D′(s) = D(s) * sin(ω(s))

	Mdiag_lo_fcn(p) = p ∈ 1:P ? h[p] * (R′(s[p]) + 2R′((s[p-1] + s[p])/2)) : 0.0
	Mdiag_hi_fcn(p) = p ∈ 0:P-1 ? h[p+1] * (R′(s[p]) + 2R′((s[p] + s[p+1])/2)) : 0.0
	Sdiag_lo_fcn(p) = p ∈ 1:P ? h[p]^2 \ quadgk(D′, s[p-1], s[p])[1] : 0.0
	Sdiag_hi_fcn(p) = p ∈ 0:P-1 ? h[p+1]^2 \ quadgk(D′, s[p], s[p+1])[1] : 0.0

	Mdiag_lo = [Mdiag_lo_fcn(p) for p ∈ 0:P]
	Mdiag_hi = [Mdiag_hi_fcn(p) for p ∈ 0:P]
	Mdiag = 3\(Mdiag_lo + Mdiag_hi)
	Moffd = 3\[h[p] * R′((s[p-1] + s[p])/2) for p ∈ 1:P]
	M = SymTridiagonal(Mdiag, Moffd)

	Sdiag_lo = [Sdiag_lo_fcn(p) for p ∈ 0:P]
	Sdiag_hi = [Sdiag_hi_fcn(p) for p ∈ 0:P]
	Sdiag = 2*(Sdiag_lo + Sdiag_hi)
	Soffd = -2 * [h[p]^2 \ quadgk(D′, s[p-1], s[p])[1] for p ∈ 1:P]
	S = SymTridiagonal(Sdiag, Soffd)

	# Solving Matrix Equation
	diffuse(U, Δtn) = OffsetArray(
		(parent(M) + Δtn * parent(S)) \ (parent(M) * parent(U)),
		Origin(0)
	)

	diffuse!(U, Δtn) = push!(U, diffuse(U[end], Δtn))

	function fem_basestepsize(U₀)
		ave(a, b) = √(a*b)

		# initial bounding guess
		Δt′ = R(0.0)^2 / D(0.0)
		Δt₊ = 1e20Δt′
		Δt₋ = 1e-20Δt′

		# loop until diffussion stepsize is nicely scaled
		U₁ = diffuse(U₀, ave(Δt₋, Δt₊))
		V₀ = sum(U₀[0:pj])
		V₁ = sum(U₁[0:pj])
		ratio = 95/100
		# counter = 0
		while !isapprox(V₁, V₀*ratio, atol = 1e-4)
			if V₁ > V₀*ratio
				Δt₋ = ave(Δt₋, Δt₊)
			elseif V₁ < V₀*ratio
				Δt₊ = ave(Δt₋, Δt₊)
			end
			U₁ = diffuse(U₀, ave(Δt₋, Δt₊))
			V₀ = sum(U₀[0:pj])
			V₁ = sum(U₁[0:pj])

			# counter += 1
			# if counter > 200
			# 	@show counter
			# 	@show Δt₋, Δt₊
			# 	@show V₁ / (V₀ * ratio)
			# end
		end
		return ave(Δt₋, Δt₊)
	end

	t = OffsetArray([0.0], Origin(0))
	U = OffsetArray([H.(sj .- s)], Origin(0))
	Δt′ = fem_basestepsize(U[0])
	tol = 1e-3
	conv = false
	n = 0

	while !conv
		n += 1
		Δtn = Δt′ * (n / 30)^2
		push!(t, t[end] + Δtn)
		diffuse!(U, Δtn)
		
		err = abs.(U[end] .- u∞) |> maximum
		
		if err < tol
			conv = true
		end
		
		if maximum(abs.(U[end] .- 0.5)) > 2
			error("Unstable.")
		end
		
		if n > 5000
			@info "Taking too long. Bailing at error $(err)."
			break
		end
	end
	
	## Return results
	return t, U, pj
end