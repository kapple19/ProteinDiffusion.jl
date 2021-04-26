Nt′ = 100

function xlim(raw::Raw)
	X = 2raw.s[raw.pj]
	Xdisp = round(100X/raw.s[end], digits = 2)
	Nt = min(Nt′, length(raw.U))

	return X, Xdisp, Nt
end

@recipe function plot(raw::Raw)
	X, Xdisp, Nt = xlim(raw)

	xlims := (0, X)
	ylims := (0, 1)

	title := raw.mode * ": Raw Output"
	xguide := "Arc Length ($Xdisp% of domain)"
	yguide := "Relative Concentration"
	legend := false

	color_palette := palette(:blues, Nt)

	raw.s, raw.U[0:Nt-1]
end

function xlim(arc::Arc)
	X = 2arc.sj
	Xdisp = round(100X/arc.smax, digits = 2)
	Nt = min(Nt′, arc.tmax)

	return X, Xdisp, Nt
end

@recipe function plot(arc::Arc)
	X, Xdisp, Nt = xlim(arc)
	
	xlims := (0, X)
	ylims := (0, 1)

	title := arc.mode * ": Interpolated Solution"
	xguide := "Arc Length ($Xdisp% of domain)"
	yguide := "Relative Concentration"
	legend := false
	
	# color_palette := palette(:blues, Nt)

	[s -> arc.u(s, t) for t ∈ LinRange(0, arc.tmax, 100)]
end