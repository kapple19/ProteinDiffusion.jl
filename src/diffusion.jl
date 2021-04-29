struct RawOutput <: PD
	mode::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

struct ArcLength <: PD
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

struct Intensity <: PD
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