Nt′ = 100

# User Recipe
# Type Recipe
# * My type to plottable data
# Plot Recipe
# Series Recipe

# Type Recipe for `Membrane`
@recipe function plot(mem::Membrane)
	ϕ = LinRange(0, 2π, 101)

	x = mem.R * cos.(ϕ)
	z = mem.R * sin.(ϕ)

	x, z
end

# @userplot FusionSlice

# mutable struct FusionSlice{FM} where FM <: FusionMode
# 	fm::FM      
# end

# @recipe function plot(fs::FusionSlice{F}) where F <: FullFusion
	
# end

# @recipe function plot(fs::FusionSlice{K}) where K <: KNRFusion

# end

# @userplot MultipleSlices

# @recipe function plot(ms::MultipleSlices)

# end

# temporal_palette(n::Int) = :blue
# temporal_palette(n::Vector{Int}) = palette(:blues, length(n))

@recipe function plot(
	raw::RawOutput;
	n = 0:min(Nt′, length(raw.U))-1)

	pal = [c for c ∈ palette(:blues, Nt′)[Nt′ .- n]]
	color_palette := pal
	# seriescolor --> pal
	label --> ""
	xmax = 2raw.s[raw.pj]
	xmaxperc = round(100xmax/raw.s[end], digits = 2)
	xlims --> (0.0, xmax)
	yguide --> "Relative Concentration"
	xguide --> "Arc Length ($xmaxperc% of span)"
	title --> raw.mode * ": Raw Data\n" * "(timesteps $(n[1]) to $(n[end]))"
	raw.s, [raw.U[n]]
end

@recipe function plot(arc::ArcLength)
	t = LinRange(0, arc.tmax, Nt′)
	label --> ""
	title --> arc.mode * ": Arc Length\n" * "($Nt′ equidistant timesteps)"
	color_palette := palette(:blues, length(t))
	xmax = 2arc.sj
	xmaxperc = round(2arc.sj/arc.smax, digits = 2)
	xlims --> (0.0, xmax)
	xguide --> "Arc Length ($xmaxperc% of span)"
	yguide --> "Relative Concentration"
	[s -> arc.u(s, t′) for t′ ∈ t]
end

@recipe function plot(int::Intensity)
	membranes = [:u, :v, :c]
	label --> ["Total" "Vesicle" "Cell"]
	title --> int.mode * ": Integrated Concentration"
	xguide --> "Time"
	yguide --> "Integrated Concentration"
	xlims --> (0, int.tmax)
	[getproperty(int, m) for m ∈ membranes]
end