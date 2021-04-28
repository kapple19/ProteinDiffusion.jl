function diffusion_fem(
	sj::Float64,
	sP::Float64,
	ω::Function,
	R::Function,
	D::Function,
	u∞::Float64)

	# Spatial Grid
	P = 1500
	pj = P ÷ 2
	ℙ = 0:P

	function spatial_grid(p::Integer)
		p ∉ ℙ && error("Index outside grid.")
		p ≤ pj && return sj * (1 - (1 - p / pj)^3)
		p > pj && return sj + (sP - sj) * ((p - pj) / (P - pj))^3
		return NaN
	end

	s = OffsetArray(
		[spatial_grid(p) for p ∈ ℙ],
		Origin(0)
	)
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
	diffuse!(U, Δtn) = push!(
		U, OffsetArray(
			(parent(M) + Δtn * parent(S)) \ (parent(M) * parent(U[end])),
			Origin(0)
		)
	)

	R²(s) = R(s)^2
	Δt′ = R²(0.0) * R²(sP) / D(0.0) / D(sP) / 1e4
	tol = 1e-3
	t = OffsetArray([0.0], Origin(0))
	U = OffsetArray([H.(sj .- s)], Origin(0))
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
		
		if maximum(U[end] .- 0.5) > 2
			error("Unstable.")
		end
		
		if n > 5000
			@info "Taking too long. Bailing at error $(err)."
			break
		end
	end
	
	## Return results
	return s, t, U, pj
end