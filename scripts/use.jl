##
using ProteinDiffusion
using Plots

##
v = Membrane(1.0, 1.0)
c = Membrane(1.0, 1.0)
Rj = 1.0

f = FullFusion(v, c)
k = KNRFusion(v, c, Rj)

##
plot(f.raw)
