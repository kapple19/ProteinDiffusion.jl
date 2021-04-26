function full_fusion(
	Rv::Float64,
	Rc::Float64,
	Dv::Float64,
	Dc::Float64)

	Rv² = Rv^2
	Rc² = Rc^2
	R′² = Rv² + Rc²
	R′ = √R′²
	Rs = 2Rv*Rc/R′
	ϕj = acos((Rc² - Rv²)/R′²)

	sj = R′ * ϕj
	sP = R′ * π
	ω(s) = s / R′
	R(s) = R′
	D(s) = Dv * H(sj - s) + Dc * H(s - sj)

	u∞ = (1 - cos(ϕj))/2

	s, t, U, pj = diffusion_fem(sj, sP, ω, R, D, u∞)

	raw = Raw("Full Fusion", s, t, U, pj)

	return Diffusion(raw.fusion, raw)
end

function knr_fusion(
	Rv::Float64,
	Rc::Float64,
	Rj::Float64,
	Dv::Float64,
	Dc::Float64)

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
	ω(s) = φ(s) * H(sj - s) + ψ(s) * H(s - sj)
	D(s) = Dv * H(sj - s) + Dc * H(s - sj)
	R(s) = Rv′ * H(sj - s) + Rc′ * H(s - sj)
	
	u∞ = Rv′² * (1 - cos(φv)) / (Rv′² * (1 - cos(φv)) + Rc′² * (1 - cos(ψc)))

	s, t, U, pj = diffusion_fem(sj, sP, ω, R, D, u∞)

	raw = Raw("KNR Fusion", s, t, U, pj)

	return Diffusion(raw.fusion, raw)
end

function fusion(
	Rv::Float64,
	Rc::Float64,
	Rj::Float64,
	Dv::Float64,
	Dc::Float64)

	f = full_fusion(Rv, Rc, Dv, Dc)
	k = knr_fusion(Rv, Rc, Rj, Dv, Dc)

	return f, k
end