function fem_diffusion(
	φ::Function,
	ψ::Function,
	D::Function,
	R::Function,
	sj::Float64,
	sP::Float64,
	uf)

	P = 1000
	pj = P÷3
	ℙ = Int.(0:P)
	s = OffsetArray(
		[
			LinRange(0.0, sj, pj+1);
			LinRange(sj, sP, P - pj + 1)
		] |> unique |> sort,
		Origin(0)
	)
	h = s |> parent |> diff

	ω(s) = φ(s) * H(sj - s) + ψ(s) * H(s - sj)
	R′(s) = R(s) * sin(ω(s))
	D′(s) = D(s) * sin(ω(s))

	Mdiag_lo_fcn(p) = p ∈ 1:P ? h[p] * (R′(s[p]) + 2R′((s[p-1] + s[p])/2)) : 0.0
	Mdiag_lo = [Mdiag_lo_fcn(p) for p ∈ 0:P]
	Mdiag_hi_fcn(p) = p ∈ 0:P-1 ? h[p+1] * (R′(s[p]) + 2R′((s[p] + s[p+1])/2)) : 0.0
	Mdiag_hi = [Mdiag_hi_fcn(p) for p ∈ 0:P]
	Mdiag = 3\(Mdiag_lo + Mdiag_hi)
	Moffd = 3\[h[p] * R′((s[p-1] + s[p])/2) for p ∈ 1:P]
	M = SymTridiagonal(Mdiag, Moffd)

	Sdiag_lo_fcn(p) = p ∈ 1:P ? h[p]^2 \ quadgk(D′, s[p-1], s[p])[1] : 0.0
	Sdiag_lo = [Sdiag_lo_fcn(p) for p ∈ 0:P]
	Sdiag_hi_fcn(p) = p ∈ 0:P-1 ? h[p+1]^2 \ quadgk(D′, s[p], s[p+1])[1] : 0.0
	Sdiag_hi = [Sdiag_hi_fcn(p) for p ∈ 0:P]
	Sdiag = 2*(Sdiag_lo + Sdiag_hi)
	Soffd = -2 * [h[p]^2 \ quadgk(D′, s[p-1], s[p])[1] for p ∈ 1:P]
	S = SymTridiagonal(Sdiag, Soffd)

	diffuse!(U, Δtn) = push!(
		U, OffsetArray(
			(parent(M) + Δtn * parent(S)) \ (parent(M) * parent(U[end])),
			Origin(0)
		)
	)

	Δt′ = R(0.0) * R(0.0) / D(0.0) / D(sP) / 1e3
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
		
		err = abs.(U[end] .- uf) |> maximum
		
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
	
	pv = 0:pj
	pc = pj:P
	vgrid = OffsetArray(s[pv], Origin(0))
	cgrid = OffsetArray(s[pc], Origin(0))
	Umat = OffsetArray(hcat(U...), Origin(0))
	Vmat = OffsetArray(Umat[pv, :], Origin(0))
	Cmat = OffsetArray(Umat[pc, :], Origin(0))

	function membrane_interpolate(grid, mat)
		itp = interpolate(
			(grid |> parent, t |> parent),
			mat |> parent,
			Gridded(Linear())
		)
		u(a, t) = grid[begin] ≤ a ≤ grid[end] ? itp(a, t) : 0.0
	end

	vs = membrane_interpolate(vgrid, Vmat)
	cs = membrane_interpolate(cgrid, Cmat)

	V = [Vmat[:, n] for n ∈ eachindex(Vmat[1, :])]
	C = [Cmat[:, n] for n ∈ eachindex(Cmat[1, :])]

	return U, V, C, s, t, vs, cs, t[end]
end