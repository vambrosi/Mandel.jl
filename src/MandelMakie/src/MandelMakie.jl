module MandelMakie

export Viewer, Viewer3D

using GLMakie, Symbolics, Parameters

import Base: show

include("Dynamics.jl")
include("Parser.jl")
include("Viewer2D.jl")
include("Viewer3D.jl")

end