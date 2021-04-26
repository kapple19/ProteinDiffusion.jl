## Load Package
using ProteinDiffusion
using Plots

cell_pars = (
	unreal = (
		Rv = 1.0,
		Rc = 2.0,
		Rj = 0.4,
		Dv = 1.0,
		Dc = 0.2
	),
	beta = (
		Rv = 150e-3, # [μm] insulin vesicles
		Rc = 4.0, # [μm] β-cells
		Rj = 50e-3,
		Dv = 1.0,
		Dc = 0.2
	),
	adipocyte = (
		Rv = 75e-3, # [μm] GLUT4 vesicles
		Rc = 17.0, # [μm] adipocytes (smaller of bimodal size distribution)
		Rj = 60e-3, # [μm] pore junction radius
		Dv = 1.0,
		Dc = 0.2
	)	
)

function run_demo(cell_type::Symbol)
	Rv, Rc, Rj, Dv, Dc = getproperty(cell_pars, cell_type)

	f, k = fusion(Rv, Rc, Rj, Dv, Dc)

	pfraw = plot(f.raw)
	pkraw = plot(k.raw)

	pfraw |> display
	pkraw |> display
end

## Run Examples
run_demo(:unreal)
run_demo(:beta)
run_demo(:adipocyte)

##
