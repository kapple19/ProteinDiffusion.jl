abstract type DiffusionSolution <: PD end
abstract type DiffusionMode <: DiffusionSolution end
abstract type DiffusionAngleSolution <: DiffusionSolution end

struct RawOutput <: DiffusionSolution
	mode::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

struct ArcLength <: DiffusionSolution
	mode::String
	u::Function
	v::Function
	c::Function
	sj::Float64
	smax::Float64
	tmax::Float64

	function ArcLength(raw::RawOutput)
		s = parent(raw.s)
		t = parent(raw.t)
		U = hcat(raw.U...)

		itp = interpolate(
			(s, t), U,
			Gridded(Linear())
		)

		sj = raw.s[raw.pj]
		smax = raw.s[end]
		tmax = raw.t[end]

		function u(s, t)
			s ∉ 0.0..smax && return 0.0
			t ∉ 0.0..tmax && return NaN
			return itp(s, t)
		end

		v(s::Real, t::Real) = 0 ≤ s ≤ sj ? u(s, t) : 0.0
		c(s::Real, t::Real) = sj ≤ s ≤ smax ? u(s, t) : 0.0

		return new(raw.mode, u, v, c, sj, smax, raw.t[end])
	end
end

struct AngleFC <: DiffusionAngleSolution
	mode::String
	u::Function
	v::Function
	c::Function
	ϕj::Float64
	tmax::Float64

	function AngleFC(fus::FusionFC, arc::ArcLength)
		ϕ2s(ϕ) = fus.R * ϕ
		u(ϕ, t) = arc.u(ϕ2s(ϕ), t)
		v(ϕ, t) = arc.v(ϕ2s(ϕ), t)
		c(ϕ, t) = arc.c(ϕ2s(ϕ), t)
		ϕj = arc.sj / fus.R
		return new(arc.mode, u, v, c, ϕj, arc.tmax)
	end
end

struct AngleKR <: DiffusionAngleSolution
	mode::String
	v::Function
	c::Function
	φmax::Float64
	ψmin::Float64
	tmax::Float64

	function AngleKR(fus::FusionKR, arc::ArcLength)
		φ2s(φ) = fus.Rv * φ
		ψ2s(ψ) = arc.sj + fus.Rc * (ψ + fus.ψc - π)
		
		v(φ, t) = arc.v(φ2s(φ), t)
		c(ψ, t) = arc.c(ψ2s(ψ), t)
		
		return new(arc.mode, v, c, fus.φv, π - fus.ψc, arc.tmax)
	end
end

struct Intensity <: DiffusionSolution
	mode::String
	u::Function
	v::Function
	c::Function
	tmax::Float64

	function Intensity(
		raw::RawOutput,
		arc::ArcLength,
		R::Function,
		ω::Function)

		pj = raw.pj
		P = lastindex(raw.s)
		N = lastindex(raw.t)

		Uint = OffsetArray(
			[
				OffsetArray(
					[
						R(raw.s[p]) * sin(ω(raw.s[p])) * raw.U[n][p]
						for p ∈ 0:P
					],
					Origin(0)
				) for n ∈ 0:N
			],
			Origin(0)
		)

		V = [integrate(raw.s[0:pj], Uint[n][0:pj]) for n ∈ eachindex(raw.t)]
		C = [integrate(raw.s[pj:P], Uint[n][pj:P]) for n ∈ eachindex(raw.t)]
		
		function itp(t, U)
			itp_ = LinearInterpolation(
				raw.t |> parent,
				U |> parent
			)

			return itp_(t)
		end

		v(t) = itp(t, V)
		c(t) = itp(t, C)
		u(t) = v(t) + c(t)

		return new(raw.mode, u, v, c, raw.t[end])
	end
end

struct DiffusionFC <: DiffusionMode
	fus::FusionFC
	raw::RawOutput
	arc::ArcLength
	ang::AngleFC
	int::Intensity

	function DiffusionFC(fus::FusionFC)
		sj = fus.R * fus.ϕj
		sP = π * fus.R

		ω(s) = s / fus.R
		R(s) = fus.R
		D(s) = fus.ves.D * H(sj - s) + fus.cel.D * H(s - sj)

		u∞ = (1 - cos(fus.ϕj))/2

		s, t, U, pj = diffusion_fem(sj, sP, ω, R, D, u∞)
		
		raw = RawOutput("Full-Collapse Fusion", s, t, U, pj)
		arc = ArcLength(raw)
		ang = AngleFC(fus, arc)
		int = Intensity(raw, arc, R, ω)

		return new(fus, raw, arc, ang, int)
	end
end

DiffusionFC(args...) = FusionFC(args...) |> DiffusionFC

struct DiffusionKR <: DiffusionMode
	fus::FusionKR
	raw::RawOutput
	arc::ArcLength
	ang::AngleKR
	int::Intensity

	function DiffusionKR(fus::FusionKR)
		sj = fus.Rv * fus.φv
		sP = sj + fus.Rc * fus.ψc
		
		φ(s) = s / fus.Rv
		ψ(s) = (s - sj) / fus.Rc + π - fus.ψc
		ω(s) = φ(s) * H(sj - s) + ψ(s) * H(s - sj)
		D(s) = fus.ves.D * H(sj - s) + fus.cel.D * H(s - sj)
		R(s) = fus.Rv * H(sj - s) + fus.Rc * H(s - sj)
		
		u∞ = fus.Rv^2 * (1 - cos(fus.φv)) / (
			fus.Rv^2 * (1 - cos(fus.φv))
			+ fus.Rc^2 * (1 - cos(fus.ψc))
		)

		s, t, U, pj = diffusion_fem(sj, sP, ω, R, D, u∞)

		raw = RawOutput("Kiss-and-Run Fusion", s, t, U, pj)
		arc = ArcLength(raw)
		ang = AngleKR(fus, arc)
		int = Intensity(raw, arc, R, ω)

		return new(fus, raw, arc, ang, int)
	end
end

DiffusionKR(args...) = FusionKR(args...) |> DiffusionKR