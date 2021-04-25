## 
using ProteinDiffusion

##
Rv = 1.0
Rc = 2.0
Rj = 0.4
Dv = 1.0
Dc = 0.2

f = full_fusion(Rv, Rc, Dv, Dc)
k = knr_fusion(Rv, Rc, Rj, Dv, Dc)

##
using Plots

anim = @animate for t ∈ LinRange(0.0, f.ang.tmax, 101)
		pf = plot(
			ϕ -> f.ang.u(ϕ, t),
			xlims = (0, π),
			ylims = (0, 1),
			title = "Full Fusion",
			ylabel = "Concentration Level",
			xlabel = "Polar Angle [rad]",
			legend = false
		)
		pk = plot(
			s -> k.arc.u(s, t),
			xlims = (0, k.arc.S),
			ylims = (0, 1),
			title = "KNR Fusion",
			xlabel = "Arc Length",
			legend = false
		)
		plot(pf, pk, layout = (1, 2))
end

gif(anim, "img/generic_2danim.gif")
