struct Membrane <: PD
	R::Float64
	D::Float64
end

struct FullFusion <: PD
	v::Membrane
	c::Membrane
	raw::RawOutput
	arc::ArcLength
	int::Intensity

	function FullFusion(v::Membrane, c::Membrane)
		Rv = v.R
		Rc = c.R
		Dv = v.D
		Dc = c.D

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

		raw = RawOutput("Full Fusion", s, t, U, pj)
		arc = ArcLength(raw)
		int = Intensity(raw, arc, R, ω)

		return new(v, c, raw, arc, int)
	end
end

struct KNRFusion <: PD
	v::Membrane
	c::Membrane
	Rj::Float64
	raw::RawOutput
	arc::ArcLength
	int::Intensity

	function KNRFusion(v::Membrane, c::Membrane, Rj::Float64)
		Rv = v.R
		Rc = c.R
		Dv = v.D
		Dc = c.D

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

		raw = RawOutput("KNR Fusion", s, t, U, pj)
		arc = ArcLength(raw)
		int = Intensity(raw, arc, R, ω)

		return new(v, c, Rj, raw, arc, int)
	end
end