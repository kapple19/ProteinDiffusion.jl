Nt′ = 100

@recipe function f(raw::Raw)
	X = 2raw.s[raw.pj]
	Xdisp = round(100X/raw.s[end], digits = 2)
	Nt = min(Nt′, length(raw.U))

	xlims := (0, X)
	ylims := (0, 1)

	title := raw.fusion * ": Raw Data (first $Nt timesteps)"
	xguide := "Arc Length ($Xdisp% of domain)"
	yguide := "Relative Concentration"
	legend := false

	color_palette := palette(:blues, Nt)
	colorbar = true

	raw.s, raw.U[0:Nt-1]
end