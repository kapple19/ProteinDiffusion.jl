abstract type DiffusionSolution <: PD end

struct RawOutput <: DiffusionSolution
	mode::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

displayname(::RawOutput) = "Raw Output"

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

displayname(::ArcLength) = "Arc-Length Solution"

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

displayname(::Intensity) = "Intensity"

struct DiffusionFC <: DiffusionSolution
	fus::FusionFC
	raw::RawOutput
	arc::ArcLength
	int::Intensity

	function DiffusionFC(fus::FusionFC)
		sj = fus.R * fus.ϕj
		sP = fus.R * π

		ω(s) = s / fus.R
		R(s) = fus.R
		D(s) = fus.ves.D * H(sj - s) + fus.cel.D * H(s - sj)

		u∞ = (1 - cos(fus.ϕj))/2

		s, t, U, pj = diffusion_fem(sj, sP, ω, R, D, u∞)

		raw = RawOutput("Full-Collapse Fusion", s, t, U, pj)
		arc = ArcLength(raw)
		int = Intensity(raw, arc, R, ω)

		return new(fus, raw, arc, int)
	end
end

DiffusionFC(args...) = FusionFC(args...) |> DiffusionFC

struct DiffusionKR <: DiffusionSolution
	fus::FusionKR
	raw::RawOutput
	arc::ArcLength
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
		int = Intensity(raw, arc, R, ω)

		return new(fus, raw, arc, int)
	end
end

DiffusionKR(args...) = FusionKR(args...) |> DiffusionKR