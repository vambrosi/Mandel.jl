module MandelMakie

export Viewer, Viewer3D, set_parameter!

using GLMakie, Symbolics, Parameters

import Base: show

include("Dynamics.jl")
include("Parser.jl")
include("Viewer2D.jl")
include("Viewer3D.jl")
include("Julia3D.jl")

end