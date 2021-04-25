NestedOSV = OffsetVector{OffsetVector{Float64, Vector{Float64}}, Vector{OffsetVector{Float64, Vector{Float64}}}}

struct FullDiffusionRaw <: PD
	U::NestedOSV
	V::NestedOSV
	C::NestedOSV
	ϕ::OffsetVector{Float64}
	t::OffsetVector{Float64}
end

struct KNRDiffusionRaw <: PD
	U::NestedOSV
	V::NestedOSV
	C::NestedOSV
	s::OffsetVector{Float64}
	t::OffsetVector{Float64}
end

struct FullDiffusionAngle <: PD
	u::Function
	v::Function
	c::Function
	ϕj::Float64
	tmax::Float64
end

struct KNRDiffusionArc <: PD
	u::Function
	v::Function
	c::Function
	sj::Float64
	S::Float64
	tmax::Float64
end

abstract type PDIntensity <: PD end

struct FullDiffusionIntensity <: PDIntensity
	I::Function
	v::Function
	c::Function
	tmax::Float64
end

struct KNRDiffusionIntensity <: PDIntensity
	I::Function
	v::Function
	c::Function
	tmax::Float64
end

struct FullDiffusion <: PD
	raw::FullDiffusionRaw
	ang::FullDiffusionAngle
	int::FullDiffusionIntensity
end

struct KNRDiffusion <: PD
	raw::KNRDiffusionRaw
	arc::KNRDiffusionArc
	int::KNRDiffusionIntensity
end

struct Diffusion <: PD
	f::FullDiffusion
	k::KNRDiffusion
end
