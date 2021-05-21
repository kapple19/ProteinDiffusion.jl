##
using ProteinDiffusion
using Plots

##
v = Membrane(1.0, 1.0)
c = Membrane(2.0, 0.2)
Rj = 0.4

fc = DiffusionFC(v, c)
# kr = DiffusionKR(v, c, Rj)

##
# plot(fc.ves) |> display
# plot(fc.cel) |> display
plot(fc.raw) |> display
plot(fc.arc) |> display
plot(fc.int) |> display

# plot(kr.ves) |> display
# plot(kr.cel) |> display
plot(kr.raw) |> display
plot(kr.arc) |> display
plot(kr.int) |> display

##
plot(fc.fus) |> display
plot(kr.fus) |> display

##
plot(kr.fus)
densityslicekrv!(kr, 2.0)
densityslicekrc!(kr, 2.0)

##
plot(fc.ewm.spot, xlims = (0, fc.ewm.tmax))

##
plot(fc.ewm.spot)
plot(fc.ewm.ring)