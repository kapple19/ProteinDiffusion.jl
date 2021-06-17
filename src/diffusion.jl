abstract type DiffusionSolution <: PD end
abstract type DiffusionMode <: DiffusionSolution end
abstract type DiffusionAngleSolution <: DiffusionSolution end

"""
Struct storing the raw output of the solution method implementation. Generic for both fusion modes.
"""
struct RawOutput <: DiffusionSolution
	mode::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

"""
Struct storing the bivariate function with respect to arc-length that interpolates the raw output solution. Generic for both fusion modes.
"""
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
			0.0 ≤ s ≤ smax && return itp(s, t)
			0.0 ≤ t ≤ tmax && return NaN
			return 0.0
		end

		v(s::Real, t::Real) = 0 ≤ s ≤ sj ? u(s, t) : 0.0
		c(s::Real, t::Real) = sj ≤ s ≤ smax ? u(s, t) : 0.0

		return new(raw.mode, u, v, c, sj, smax, raw.t[end])
	end
end

"""
Struct storing the Full-Collapse bivariate function interpolated from the raw output and transformed from arc-length to polar angle.
"""
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

"""
Struct storing the Kiss-and-Run bivariate function interpolated from the raw output and transformed from arc-length to polar angle.
"""
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

"""
Struct for storing the time-varying intensities. Generic for both fusion modes.

The intensity is obtained via integration of the raw solution over space for respective regions.
"""
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

# decay_fcn_cut(ζ, δ) = H(δ - ζ)
# decay_fcn_exp(ζ, δ) = exp(-ζ/δ)

# struct Decay
# 	spot::Function
# 	ring::Function

# 	function Decay(X::Float64, x2ζ⁻::Function, x2ζ⁺::Function)
# 		spot_x_bnds = (0, X)
# 		ring_x_bnds = (X, 2X)

# 		spot = get_decay("cut", spot_x_bnds, x2ζ⁻, x2ζ⁺)
# 		ring = get_decay("cut", ring_x_bnds, x2ζ⁻, x2ζ⁺)

# 		return new(spot, ring)
# 	end
# end

# function get_intensity(
# 	decay::Function,
# 	bnds::NTuple{2, Float64}, x, t, U,
# 	x2ζ::Function)

# 	@assert x |> issorted
# 	@assert length(x) == length(unique(x))

# 	x_lb = 0.0
# 	x_ub = fus.ves.R

# 	if x_lb ∉ x
# 		Uitp = [LinearInterpolate(x, U′, extrapolation_bc = 0.0) for U′ ∈ U]
# 		U_lb = [Uitp′(x_lb) for Uitp′ ∈ Uitp]
# 		for (n, U′) ∈ enumerate(U)
# 			pushfirst!(U′, U_lb[n])
# 		end
# 		push!(x, x_lb)
# 		sort!(x)
# 	end

# 	if x_ub ∉ x
# 		Uitp = [LinearInterpolate(x, U′, extrapolation_bc = 0.0) for U′ ∈ U]
# 		U_ub = [Uitp′(x_ub) for Uitp′ ∈ Uitp]
# 		for (n, U′) ∈ enumerate(U)
# 			push!(U′, U_ub[n])
# 		end
# 		push!(x, x_ub)
# 		sort!(x)
# 	end

# 	p = find(x_lb .≤ x .≤ x_ub)

# 	x = x[p]
# 	U = [U′[p] for U′ ∈ U]

# 	Imat = [integrate(x, U′) * decay(x2ζ(x), δ) for U′ ∈ U]

# 	Iitp = LinearInterpolation(t, Imat)

# 	I(t) = Iitp(t)

# 	return I
# end

# function get_decay(
# 	decay::Function,
# 	x_bnds::NTuple{2, Float64},
# 	x2ζ⁻::Function,
# 	x2ζ⁺::Function)

# 	I⁻ = get_intensity(decay, x_bnds, x⁻, t, U⁻, x2ζ⁻)
# 	I⁺ = get_intensity(decay, x_bnds, x⁺, t, U⁺, x2ζ⁺)

# 	I(t) = I⁻(t) + I⁺(t)
# end

# get_decay(decay::String, args...) = get_decay(
# 	getfield(ProteinDiffusion, "decay_fcn_" * decay),
# 	args...
# )

# """
# Struct for storing the evanescent wave microscopy intensity model functions obtained from transformation of the raw data to viewing range and depth and integrating over pre-defined viewing range range and factoring each concentration level with respect to depth.

# Separated from the `Intensity` struct due to difference inherent difference in integration domain.
# """
# struct Microscopy <: DiffusionSolution
# 	mode::String
# 	cut::Decay
# 	exp::Decay
# 	xmax::Float64
# 	tmax::Float64
# end

# function Microscopy(fus::FusionFC, raw::RawOutput)
# 	decay_cut = Decay("cut", fus.R)
# 	decay_exp = Decay("exp", fus.R)

# 	Microscopy(get_fusion_mode(fus), decay_cut, decay_exp, raw.t[end])
# end

# function Microscopy(fus::FusionKR, raw::RawOutput)
	

# 	Microscopy(get_fusion_mode(fus), decay_cut, decay_exp, raw.t[end])
# end

"""
Struct for storing the diffusion information for Full-Collapse Fusion, including evanescent wave microscopy.
"""
struct DiffusionFC <: DiffusionMode
	fus::FusionFC
	raw::RawOutput
	arc::ArcLength
	ang::AngleFC
	int::Intensity
	# ewm::Microscopy

	function DiffusionFC(
		fus::FusionFC,
		δ::Float64 = Inf64)

		sj = fus.R * fus.ϕj
		sP = π * fus.R

		ω(s) = s / fus.R
		R(s) = fus.R
		D(s) = fus.ves.D * H(sj - s) + fus.cel.D * H(s - sj)

		u∞ = (1 - cos(fus.ϕj))/2

		P′ = 1500
		pj′ = P′ ÷ 2

		function spatial_grid(p::Int64)
			p ∉ 0:P′ && error("Index outside grid.")
			p == 0 && return 0.0
			p == pj′ && return sj
			p == P′ && return sP
			p < pj′ && return sj * (1 - (1 - p / pj′)^3)
			p > pj′ && return sj + (sP - sj) * ((p - pj′) / (P′ - pj′))^3
			return NaN
		end

		s = OffsetVector(
			[
				[spatial_grid(p) for p ∈ 0:P′];
				π/4 * fus.R
			] |> unique! |> sort!,
			Origin(0)
		)
		pj = findfirst(s .== sj)

		t, U, pj = diffusion_fem(s, pj, ω, R, D, u∞)
		
		raw = RawOutput(get_fusion_mode(fus), s, t, U, pj)
		arc = ArcLength(raw)
		ang = AngleFC(fus, arc)
		int = Intensity(raw, arc, R, ω)
		# ewm = Microscopy(fus, raw, δ)

		# return new(fus, raw, arc, ang, int, ewm)
		return new(fus, raw, arc, ang, int)
	end
end

DiffusionFC(args...) = FusionFC(args...) |> DiffusionFC

"""
Struct for storing the diffusion information for Kiss-and-Run Fusion, including evanescent wave microscopy.
"""
struct DiffusionKR <: DiffusionMode
	fus::FusionKR
	raw::RawOutput
	arc::ArcLength
	ang::AngleKR
	int::Intensity
	# ewm::Microscopy

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

		P′ = 1500
		pj′ = P′ ÷ 2

		φ2s(φ) = φ * fus.Rv
		ψ2s(ψ) = fus.Rc * (ψ + fus.ψc - π) + sj

		function spatial_grid(p::Int64)
			p ∉ 0:P′ && error("Index outside grid.")
			p == 0 && return 0.0
			p == pj′ && return sj
			p == P′ && return sP
			p < pj′ && return sj * (1 - (1 - p / pj′)^3)
			p > pj′ && return sj + (sP - sj) * ((p - pj′) / (P′ - pj′))^3
			return NaN
		end

		s = OffsetVector(
			[
				[spatial_grid(p) for p ∈ 0:P′];
				φ2s(π/4); ψ2s(π/4) 
			] |> unique! |> sort!,
			Origin(0)
		)
		pj = findfirst(s .== sj)

		t, U, pj = diffusion_fem(s, pj, ω, R, D, u∞)

		raw = RawOutput(get_fusion_mode(fus), s, t, U, pj)
		arc = ArcLength(raw)
		ang = AngleKR(fus, arc)
		int = Intensity(raw, arc, R, ω)
		# ewm = Microscopy(fus, raw, δ)

		# return new(fus, raw, arc, ang, int, ewm)
		return new(fus, raw, arc, ang, int)
	end
end

DiffusionKR(args...) = FusionKR(args...) |> DiffusionKR

get_fusion_mode(::FusionFC) = "Full-Collapse"
get_fusion_mode(::FusionKR) = "Kiss-and-Run"
get_fusion_mode(dm::DiffusionMode) = get_fusion_mode(dm.fus)
get_fusion_mode(ds::DiffusionSolution) = ds.mode