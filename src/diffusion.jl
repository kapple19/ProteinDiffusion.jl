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

# TODO: Replace quadgk on interpolated fcns with integrate on the raw data

struct Microscopy <: DiffusionSolution
	U::Function
	ring::Function
	spot::Function
	xmax::Float64
	tmax::Float64

	function Microscopy(U::Function, Rv::Float64, xmax::Float64, tmax::Float64)
		Uring(t) = quadgk(x -> U(x, t), 0.0, Rv)[1]
		Uspot(t) = quadgk(x -> U(x, t), Rv, 2Rv)[1]

		return new(U, Uring, Uspot, xmax, tmax)
	end
end

function Microscopy(δ::Float64, fus::FusionFC, ang::AngleFC)
	function ϕ(x)
		ϕ′ = asin(x / fus.R)
		return sort([ϕ′ π-ϕ′], dims = 2)
	end

	z(ϕ) = fus.R * cos(ϕ)

	Usummagrand(x, t) = ang.u.(ϕ(x), t) .* exp.(-z.(ϕ(x)) / δ)
	U(x, t) = sum(Usummagrand(x, t), dims = 2)[1]

	return Microscopy(U, fus.ves.R, fus.R, ang.tmax)
end

function Microscopy(δ::Float64, fus::FusionKR, ang::AngleKR)
	function φ(x)
		φ′ = asin(x / fus.Rv)
		φset = Float64[]
		0.0 ≤ φ′ ≤ fus.φv && push!(φset, φ′)
		0.0 ≤ π - φ′ ≤ fus.φv && push!(φset, π - φ′)
		return sort(φset', dims = 2)
	end

	function ψ(x)
		ψ′ = asin(x / fus.Rc)
		ψset = Float64[]
		π - fus.ψc ≤ ψ′ ≤ π && push!(ψset, ψ′)
		π - fus.ψc ≤ π - ψ′ ≤ π && push!(ψset, ψ′)
		return sort(ψset', dims = 2)
	end

	ξ(φ) = fus.Rv * cos(φ)
	ζ(ψ) = fus.Rc * cos(ψ)

	Vsummagrand(x, t) = ang.v.(φ(x), t) .* exp.(-ξ.(φ(x)) / δ)
	Csummagrand(x, t) = ang.c.(ψ(x), t) .* exp.(-ζ.(ψ(x)) / δ)

	V(x, t) = sum(Vsummagrand(x, t), dims = 2)[1]
	C(x, t) = sum(Csummagrand(x, t), dims = 2)[1]
	U(x, t) = V(x, t) + C(x, t)

	return Microscopy(U, fus.ves.R, fus.Rc, ang.tmax)
end

struct DiffusionFC <: DiffusionMode
	fus::FusionFC
	raw::RawOutput
	arc::ArcLength
	ang::AngleFC
	int::Intensity
	ewm::Microscopy

	function DiffusionFC(
		fus::FusionFC,
		δ::Float64 = Inf64
	)
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
		ewm = Microscopy(δ, fus, ang)

		return new(fus, raw, arc, ang, int, ewm)
	end
end

DiffusionFC(args...) = FusionFC(args...) |> DiffusionFC

struct DiffusionKR <: DiffusionMode
	fus::FusionKR
	raw::RawOutput
	arc::ArcLength
	ang::AngleKR
	int::Intensity
	ewm::Microscopy

	function DiffusionKR(
		fus::FusionKR,
		δ::Float64 = Inf
	)
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
		ewm = Microscopy(δ, fus, ang)

		return new(fus, raw, arc, ang, int, ewm)
	end
end

DiffusionKR(args...) = FusionKR(args...) |> DiffusionKR