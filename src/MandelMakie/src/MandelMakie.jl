"""
MandelMakie.jl is a module to explore complex dynamical systems. Its main functions are:

- `Viewer(f)`: plots the Mandelbrot and Julia (or just Julia) sets associated with `f`.
- `Viewer3D(f)`: plots the same sets but on the 2-sphere (`f` must be rational).
"""
module MandelMakie

export Viewer, Viewer3D, set_parameter!

using GLMakie, Symbolics, Parameters

import Base: show

include("Dynamics.jl")
include("Parser.jl")
include("Viewer2D.jl")
include("Viewer3D.jl")
include("Fatou3D.jl")

end
