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
		Rv = 75.0, # [μm] GLUT4 vesicles
		Rc = 17e3, # [μm] adipocytes (smaller of bimodal size distribution)
		Rj = 60.0, # [μm] pore junction radius
		Dv = 1.0,
		Dc = 0.2
	)	
)

function run_demo(CellType::Symbol)
	Rv, Rc, Rj, Dv, Dc = getproperty(cell_pars, CellType)

	v = Membrane(Rv, Dv)
	c = Membrane(Rc, Dc)

	f = FullFusion(v, c)
	k = KNRFusion(v, c, Rj)

	pfraw = plot(f.raw, title = "Full Fusion: Raw Data\n" * string(CellType) * " Cell")
	pfarc = plot(f.arc, title = "Full Fusion: Arc Length\n" * string(CellType) * " Cell")
	pfint = plot(f.int, title = "Full Fusion: Integrated Concentration\n" * string(CellType) * " Cell")

	pkraw = plot(k.raw, title = "KNR Fusion: Raw Data\n" * string(CellType) * " Cell")
	pkarc = plot(k.arc, title = "KNR Fusion: Arc Length\n" * string(CellType) * " Cell")
	pkint = plot(k.int, title = "KNR Fusion: Integrated Concentration\n" * string(CellType) * " Cell")

	display.(
		[
			pfraw,
			pfarc,
			pfint,
			pkraw,
			pkarc,
			pkint
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

# Run Examples
run_demo(:Unrealistic)
run_demo(:Beta)
run_demo(:Adipocyte)
