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
R = √(Rv^2 + Rc^2)
ztoϕ(z) = acos(z/R)
uz(z, t) = f.ang.u(z |> ztoϕ, t)
z = R*LinRange(1, -1, 101)

using GLMakie

sphere = Sphere(Point3(0.0), R)

fig, scn, m = mesh(
	sphere,
	colormap = :blues,
	color = hcat(repeat([uz.(z, 0.0)], 2)...)
)

function col(t)
	c = hcat(repeat([uz.(z, t)], 2)...)
	c[end, :] = [1 1]
	c
end

T = f.ang.tmax/10
record(
	fig,
	"img/generic_3danim.gif",
	LinRange(0, T, 21),
	compression = 51,
	framerate = 6) do t
	m.color = col(t)
end
