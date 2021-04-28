## Load Package
using ProteinDiffusion
using Plots

cell_pars = (
	Unrealistic = (
		Rv = 1.0,
		Rc = 2.0,
		Rj = 0.4,
		Dv = 1.0,
		Dc = 0.2
	),
	Beta = (
		Rv = 150.0, # [μm] insulin vesicles
		Rc = 4.0e3, # [μm] β-cells
		Rj = 50.0,
		Dv = 1e3,
		Dc = 2e2
	),
	Adipocyte = (
		Rv = 75e-3, # [μm] GLUT4 vesicles
		Rc = 17.0, # [μm] adipocytes (smaller of bimodal size distribution)
		Rj = 60e-3, # [μm] pore junction radius
		Dv = 1.0,
		Dc = 0.2
	)	
)

function run_demo(CellType::Symbol)
	Rv, Rc, Rj, Dv, Dc = getproperty(cell_pars, CellType)

	f, k = fusion(Rv, Rc, Rj, Dv, Dc)

	pfraw = plot(f.raw, title = string(CellType) * " Cell")
	pfarc = plot(f.arc, title = string(CellType) * " Cell")
	@time pfint = plot(f.int, title = string(CellType) * " Cell")

	pkraw = plot(k.raw, title = string(CellType) * " Cell")
	pkarc = plot(k.arc, title = string(CellType) * " Cell")
	@time pkint = plot(k.int, title = string(CellType) * " Cell")

	display.(
		[
			pfraw, pfarc, pfint,
			pkraw, pkarc, pkint
		]
	)

	celltype = CellType |> string |> lowercase
	
	savefig(pfraw, "plots/$(celltype)_fullfusion_raw.png")
	savefig(pfarc, "plots/$(celltype)_fullfusion_arc.png")
	savefig(pfint, "plots/$(celltype)_fullfusion_int.png")
	savefig(pkraw, "plots/$(celltype)_knrfusion_raw.png")
	savefig(pkarc, "plots/$(celltype)_knrfusion_arc.png")
	savefig(pkint, "plots/$(celltype)_knrfusion_int.png")
end

## Run Examples
run_demo(:Unrealistic)
run_demo(:Beta)
run_demo(:Adipocyte)

##
Rv = 150.0
Rc = 4e3
Rj = 50.0
Dv = 1e3
Dc = 2e2

@time f = full_fusion(Rv, Rc, Dv, Dc)

##
@time plot(f.int)

##
Rv = 10.0
Rc = 10.0
Rj = 10.0
Dv = 10.0
Dc = 10.0

@time k = knr_fusion(Rv, Rc, Rj, Dv, Dc)

##
plot(k.int)

##
