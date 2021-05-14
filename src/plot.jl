Nt′ = 100

@recipe function cross_section(fc::FusionFC)
	ϕv = LinRange(-fc.ϕj, fc.ϕj, 101)
	ϕc = LinRange(fc.ϕj, 2π - fc.ϕj, 101)

	xv = @. fc.R * cos(π/2 - ϕv)
	zv = @. fc.R * sin(π/2 - ϕv)
	xc = @. fc.R * cos(π/2 - ϕc)
	zc = @. fc.R * sin(π/2 - ϕc)

	label := ["Vesicle" "Cell"]
	aspect_ratio := 1
	
	[xv, xc], [zv, zc]
end

@recipe function cross_section(kr::FusionKR)
	φv = LinRange(π - kr.φv, π + kr.φv, 101)
	ψc = LinRange(π - kr.ψc, π + kr.ψc, 101)

	xv = @. kr.Rv * cos(π/2 - φv)
	zv = @. kr.Rc * cos(π - kr.ψc) - kr.Rv * cos(π - kr.φv)	+ kr.Rv * sin(π/2 - φv)
	xc = @. kr.Rc * cos(π/2 - ψc)
	zc = @. kr.Rc * sin(π/2 - ψc)

	label := ["Vesicle" "Cell"]
	aspect_ratio := 1

	[xv, xc], [zv, zc]
end

@recipe function plot(
	raw::RawOutput;
	n = 0:min(Nt′, length(raw.U))-1)

	pal = [c for c ∈ palette(:blues, Nt′)[Nt′ .- n]] |> reverse
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