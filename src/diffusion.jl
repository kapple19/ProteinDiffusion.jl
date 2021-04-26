struct Raw <: PD
	mode::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

struct Arc <: PD
	mode::String
	u::Function
	v::Function
	c::Function
	sj::Float64
	smax::Float64
	tmax::Float64

	function Arc(raw::Raw)
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

		v(s, t) = 0 ≤ s ≤ sj ? u(s, t) : 0.0
		c(s, t) = sj ≤ s ≤ smax ? u(s, t) : 0.0

		return new(raw.mode, u, v, c, sj, smax, raw.t[end])
	end
end

struct Diffusion <: PD
	mode::String
	raw::Raw
	arc::Arc

	Diffusion(raw) = new(raw.mode, raw, Arc(raw))
end