##
using ProteinDiffusion
using Plots

##
v = Membrane(1.0, 1.0)
c = Membrane(2.0, 0.2)
Rj = 0.4

fc = DiffusionFC(v, c)
kr = DiffusionKR(v, c, Rj)

##
densityslicefc(fc, 0.0)
plot!(fc.fus)

##
anim = @animate for t ∈ LinRange(0, fc.ang.tmax, 51)
	densityslicefc(fc, t)
	plot!(fc.fus)
end
gif(anim, "anim/fc.gif", fps = 15)

##
densityslicekrv(kr, 0.0)
densityslicekrc!(kr, 0.0)
plot!(kr.fus)

##
anim = @animate for t ∈ [
	repeat([0.0], 20);
	LinRange(0, kr.ang.tmax, 51)
]
	densityslicekrv(kr, t)
	densityslicekrc!(kr, t)
	plot!(kr.fus)
end
gif(anim, "anim/kr.gif", fps = 15)
