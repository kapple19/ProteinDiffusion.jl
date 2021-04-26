## 
using ProteinDiffusion
using Plots
using RecipesBase

##
Rv = 1.0
Rc = 2.0
Rj = 0.4
Dv = 1.0
Dc = 0.2

f, k = fusion(Rv, Rc, Dv, Dc)

##
tmax = min(f.arc.tmax, k.arc.tmax)

anim = @animate for t ∈ [repeat([0.0], 31); LinRange(0.0, tmax, 101)]
	pf = plot(
		s -> f.arc.u(s, t),
		title = "Full Fusion",
		xlabel = "Arc Length",
		ylabel = "Concentration",
		xlims = (0, f.arc.smax),
		ylims = (0, 1),
		legend = false
	)

	pk = plot(
		s -> k.arc.u(s, t),
		title = "KNR Fusion",
		xlabel = "Arc Length",
		xlims = (0, k.arc.smax),
		ylims = (0, 1),
		legend = false
	)

	plot(pf, pk, layout = (1, 2))
end

gif(anim, "anim/unrealistic_2danim.gif")
