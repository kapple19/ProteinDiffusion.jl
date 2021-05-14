##
using ProteinDiffusion
using Plots

##
v = Membrane(1.0, 1.0)
c = Membrane(2.0, 0.2)
Rj = 0.4

f = FullFusion(v, c)
k = KNRFusion(v, c, Rj)

##
plot(f.v) |> display
plot(f.c) |> display
plot(f.raw) |> display
plot(f.arc) |> display
plot(f.int) |> display

plot(k.v) |> display
plot(k.c) |> display
plot(k.raw) |> display
plot(k.arc) |> display
plot(k.int) |> display
