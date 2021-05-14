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
anim = @animate for t ∈ LinRange(0, fc.ang.tmax, 501)
	densityslicefc(fc, t)
	plot!(fc.fus)
end
gif(anim, "anim/fc.gif", fps = 15)
