## Load Package
using ProteinDiffusion
using Plots

function run_demo(Rv, Rc, Rj, Dv, Dc, name::String)
	# Run Models
	pd = fusion(Rv, Rc, Rj, Dv, Dc)

	# Plot Results
	pfang = plot(pd.f.ang)
	pfint = plot(pd.f.int)
	pkarc = plot(pd.k.arc)
	pkint = plot(pd.k.int)

	display(pfang)
	display(pfint)
	display(pkarc)
	display(pkint)

	savefig(pfang, "img/" * name * "_fullfusion_ang.png")
	savefig(pfint, "img/" * name * "_fullfusion_int.png")
	savefig(pkarc, "img/" * name * "_knrfusion_arc.png")
	savefig(pkint, "img/" * name * "_knrfusion_int.png")
end

## Generic
Rv = 1.0
Rc = 2.0
Rj = 0.4
Dv = 1.0
Dc = 0.2

run_demo(Rv, Rc, Rj, Dv, Dc, "generic")

## β-cells
Rv = 150e-3 # [μm] insulin vesicles
Rc = 4.0 # [μm] β-cells
Rj = 50e-3
Dv = 1.0
Dc = 0.2

run_demo(Rv, Rc, Rj, Dv, Dc, "betacells")

## Adipocytes
Rv = 75e-3 # [μm] GLUT4 vesicles
Rc = 17.0 # [μm] adipocytes (smaller of bimodal size distribution)
Rj = Rv/2 # [μm] pore junction radius
Dv = 1.0
Dc = 0.2

run_demo(Rv, Rc, Rj, Dv, Dc, "adipocytes")
