function full_fusion(Rv, Rc, Dv, Dc)
	Rv² = Rv^2
	Rc² = Rc^2
	R′² = Rv² + Rc²
	R′ = √R′²
	Rs = 2Rv*Rc/R′
	ϕj = acos((Rc² - Rv²)/R′²)
	
	D(ϕ) = Dv * H(ϕj - ϕ) + Dc * H(ϕ - ϕj)
	
	uf = (1 - cos(ϕj))/2
	
	U, V, C, ϕ, t, vϕ, cϕ, tmax = fem_diffusion(
		ϕ -> ϕ,
		ϕ -> ϕ,
		D,
		ϕ -> R′²,
		ϕj,
		π |> Float64,
		uf
	)

	function uϕ(ϕ::Float64, t)
		ϕ < 0.0 && return NaN
		ϕ > π && return NaN
		ϕ < ϕj && return vϕ(ϕ, t)
		return cϕ(ϕ, t)
	end

	# xtoϕ(x) = x/R′ |> asin

	# ϕspot = xtoϕ(Rs)
	# ϕs = [ϕ[ϕ .≤ ϕspot]; ϕspot] |> unique |> sort
	# ϕr = [ϕ[ϕspot*3/4 .≤ ϕ .≤ ϕspot]; ϕspot] |> unique |> sort
	ϕv = [ϕ[ϕ .≤ ϕj]; ϕj] |> unique |> sort
	ϕc = [ϕ[ϕj .≤ ϕ]; ϕj] |> unique |> sort

	I(t) = R′² * quadgk(ϕ -> sin(ϕ) * uϕ(ϕ, t), ϕ...)[1]
	# Is(t) = R′² * quadgk(ϕ -> sin(ϕ) * uϕ(ϕ, t), ϕs...)[1]
	# Ir(t) = R′² * quadgk(ϕ -> sin(ϕ) * uϕ(ϕ, t), ϕr...)[1]
	Iv(t) = R′² * quadgk(ϕ -> sin(ϕ) * vϕ(ϕ, t), ϕv...)[1]
	Ic(t) = R′² * quadgk(ϕ -> sin(ϕ) * cϕ(ϕ, t), ϕc...)[1]

	raw = FullDiffusionRaw(U, V, C, ϕ, t)
	ang = FullDiffusionAngle(uϕ, vϕ, cϕ, ϕj, tmax)
	# int = FullDiffusionIntensity(I, Is, Ir, Iv, Ic, tmax)
	int = FullDiffusionIntensity(I, Iv, Ic, tmax)

	return FullDiffusion(raw, ang, int)
end

function knr_fusion(Rv, Rc, Rj, Dv, Dc)
	Rv² = Rv^2
	Rc² = Rc^2
	Rj² = Rj^2
	Rv′² = 4Rv²^2/(4Rv² - Rj²)
	Rc′² = 4Rc²^2/(4Rc² - Rj²)
	Rv′ = √Rv′²
	Rc′ = √Rc′²
	φv = π - asin(Rj/Rv′)
	ψc = π - asin(Rj/Rc′)
	sj = Rv′ * φv
	sP = sj + Rc′ * ψc
	
	φ(s) = s/Rv′
	ψ(s) = (s - sj)/Rc′ + π - ψc
	D(s) = Dv * H(sj - s) + Dc * H(s - sj)
	R(s) = Rv′ * H(sj - s) + Rc′ * H(s - sj)
	
	uf = Rv′² * (1 - cos(φv)) / (Rv′² * (1 - cos(φv)) + Rc′² * (1 - cos(ψc)))
	
	U, V, C, s, t, vs, cs, tmax = fem_diffusion(φ, ψ, D, R, sj, sP, uf)

	function us(s::Float64, t)
		s < 0.0 && return NaN
		s > sP && return NaN
		s < sj && return vs(s, t)
		return cs(s, t)
	end

	# xtoψ(x) = x/Rc′ |> asin
	# ψtos(ψ) = sj + Rc′ * (ψ - π + ψc)
	# xtos(x) = x |> xtoψ |> ψtos
	ω(s) = φ(s) * H(sj - s) + ψ(s) * H(s - sj)

	# ss = [s[s .≤ xtos(Rv′)]; xtos(Rv′)] |> unique |> sort
	sv = [s[s .≤ sj]; sj] |> unique |> sort
	sc = [s[sj .≤ s]; sj] |> unique |> sort

	I(t) = quadgk(s -> R(s) * sin(ω(s)) * us(s, t), s...)[1]
	# Is(t) = quadgk(s -> R(s) * sin(ω(s)) * us(s, t), ss...)[1]
	# Ir(t) = quadgk(s -> R(s) * sin(ω(s)) * us(s, t), sr...)[1]
	Iv(t) = quadgk(s -> R(s) * sin(ω(s)) * vs(s, t), sv...)[1]
	Ic(t) = quadgk(s -> R(s) * sin(ω(s)) * cs(s, t), sc...)[1]

	raw = KNRDiffusionRaw(U, V, C, s, t)
	arc = KNRDiffusionArc(us, vs, cs, sj, sP, tmax)
	# int = KNRDiffusionIntensity(I, Is, Ir, Iv, Ic, tmax)
	int = KNRDiffusionIntensity(I, Iv, Ic, tmax)
	
	return KNRDiffusion(raw, arc, int)
end

function fusion(Rv, Rc, Rj, Dv, Dc)
	f = full_fusion(Rv, Rc, Dv, Dc)
	k = knr_fusion(Rv, Rc, Rj, Dv, Dc)

	return Diffusion(f, k)
end
