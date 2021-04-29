@recipe function plot(mem::Membrane)
	ϕ = LinRange(0, 2π, 101)

	x = mem.R * cos.(ϕ)
	z = mem.R * sin.(ϕ)

	x, z
end
Nt′ = 100

function xlim(raw::RawOutput)
	X = 2raw.s[raw.pj]
	Xdisp = round(100X/raw.s[end], digits = 2)

	return X, Xdisp
end

@recipe function plot(
	raw::RawOutput;
	title = "",
	n = 0:min(Nt′, length(raw.U))-1)

	markerstrokewidth := 0
	markersize := 3

	X, Xdisp = xlim(raw)
	Nt = length(n)

	xlims := (0, X)
	ylims := (0, 1)

	title := raw.mode * ": Raw Output (first $Nt timesteps)" * "\n$title Cell"^sign(length(title))
	xguide := "Arc Length ($Xdisp% of domain)"
	yguide := "Relative Concentration"
	legend := false

	color_palette := palette(:blues, Nt)

	raw.s, raw.U[n]
end

function xlim(arc::ArcLength)
	X = 2arc.sj
	Xdisp = round(100X/arc.smax, digits = 2)
	Nt = min(Nt′, arc.tmax)

	return X, Xdisp, Nt
end

@recipe function plot(
	arc::ArcLength;
	title = "",
	t = LinRange(0, arc.tmax, 100))

	X, Xdisp, Nt = xlim(arc)
	
	xlims := (0, X)
	ylims := (0, 1)

	title := arc.mode * ": Interpolated Solution" * "\n$title Cell"^sign(length(title))
	xguide := "Arc Length ($Xdisp% of domain)"
	yguide := "Relative Concentration"
	legend := false
	
	color_palette := palette(:blues, length(t))

	[s -> arc.u(s, t′) for t′ ∈ t]
end

@recipe function plot(int::Intensity; title = "", membranes = [:u, :v, :c])
	xlims := (0, int.tmax)

	title := int.mode * ": Integrated Concentration" * "\n$title Cell"^sign(length(title))
	xguide := "Time"
	yguide := "Intensity"
	label := ["Total" "Vesicle" "Cell"]

	[getproperty(int, m) for m ∈ membranes]
end