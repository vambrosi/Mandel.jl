import Pkg
Pkg.instantiate()
f(z, c) = z^2 + c
crit(c) = c
include("./src/MandelMakie.jl")
MandelMakie.Viewer(f; crit = crit, mandel_diameter = 1.0)
