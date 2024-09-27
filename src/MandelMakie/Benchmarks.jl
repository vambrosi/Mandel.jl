using BenchmarkTools
using MandelMakie

pick_parameter! = MandelMakie.Refactor.pick_parameter!

function test(c, p)
    f(z, c) = z^2 + c
    viewer =
        Viewer(f, mandel_center = -0.5, convergence_criterion = c, projective_metric = p)

    julia = viewer.julia
    mandel = viewer.mandel
    d_system = viewer.d_system
    options = viewer.options

    julia.center = -0.2395058775076436 + 0.5028211721392088im
    julia.diameter = 3.5101352053743486e-5
    parameter = -0.08424621236852714 + 0.7400185007617068im

    pick_parameter!(julia, mandel, d_system, options, parameter)

    @btime MandelMakie.Refactor.update_view!($julia, $d_system, $options)
end

function tests()
    println(rpad("Criterion", 20), rpad("Metric", 14), "Benchmark")
    for c in [:escape_time, :almost_periodic, :near_attractor]
        for p in [false, true]
            print(rpad(c, 20), rpad(p ? "projective" : "plane", 12))
            test(c, p)
        end
    end
end

tests()

function profile(c, p)
    f(z, c) = z^2 + c
    viewer =
        Viewer(f, mandel_center = -0.5, convergence_criterion = c, projective_metric = p)

    julia = viewer.julia
    mandel = viewer.mandel
    d_system = viewer.d_system
    options = viewer.options

    julia.center = -0.2395058775076436 + 0.5028211721392088im
    julia.diameter = 3.5101352053743486e-5
    parameter = -0.08424621236852714 + 0.7400185007617068im

    pick_parameter!(julia, mandel, d_system, options, parameter)

    @profview MandelMakie.Refactor.update_view!(julia, d_system, options)
end

profile(:near_attractor, false)

# Julia Version 1.10.4
# Commit 48d4fd48430 (2024-06-04 10:41 UTC)
# Build Info:
#   Official https://julialang.org/ release
# Platform Info:
#   OS: macOS (arm64-apple-darwin22.4.0)
#   CPU: 8 × Apple M2
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, apple-m1)
# Threads: 4 default, 0 interactive, 2 GC (on 4 virtual cores)
# Environment:
#   JULIA_EDITOR = code

# Project MandelMakie v0.1.0
# Status `~/GitHub/Mandel/src/MandelMakie/Project.toml`
#   [a8cc5b0e] Crayons v4.1.1
#   [e9467ef8] GLMakie v0.10.11
#   [5c1252a2] GeometryBasics v0.4.11
#   [ac1192a8] HypertextLiteral v0.9.5
# ⌃ [2edaba10] Nemo v0.46.2
#   [d96e819e] Parameters v0.12.3
#   [c3e4b0f8] Pluto v0.19.46
#   [f27b6e38] Polynomials v4.0.11
#   [1e83bf80] StaticArraysCore v1.4.3
#   [0c5d862f] Symbolics v6.13.0
#   [ade2ca70] Dates
#   [37e2e46d] LinearAlgebra
# Info Packages marked with ⌃ have new versions available and may be upgradable.

# RGBA vs RGB benchmark

# PARAMETERS
# f(z, c) = z^2 + c
# viewer = Viewer(f, mandel_center = -0.5, convergence_criterion = :near_attractor)
# viewer.julia.center = -0.2395058775076436 + 0.5028211721392088im
# viewer.julia.diameter = 3.5101352053743486e-5
# viewer.julia.parameter = -0.08424621236852714 + 0.7400185007617068im

# RESULTS
# RGBA{Float64} Render time = 1.267 s (96582828 allocations: 6.28 GiB)
# RGB{Float64}  Render time = 1.262 s (96582834 allocations: 6.28 GiB)
# RGBA{Float32} Render time = 1.260 s (96583040 allocations: 6.23 GiB)
# RGB{Float32}  Render time = 1.305 s (96582828 allocations: 6.26 GiB)

# Comparisons RGBA{Float64}
# Criterion           Metric        Benchmark
# escape_time         plane         194.617 ms (5852401 allocations: 211.92 MiB)
# escape_time         projective    669.913 ms (6600070 allocations: 284.40 MiB)
# almost_periodic     plane         210.769 ms (11530615 allocations: 469.35 MiB)
# almost_periodic     projective    1.086 s (11530616 allocations: 567.03 MiB)
# near_attractor      plane         1.392 s (97918506 allocations: 6.38 GiB)
# near_attractor      projective    1.509 s (96319459 allocations: 6.29 GiB)
# Next to last option is unreasonably slower.

# After Removing type annotations (MaybeCloseBy, MaybeOrbitData)
# Criterion           Metric        Benchmark
# escape_time         plane         192.778 ms (5852401 allocations: 211.92 MiB)
# escape_time         projective    665.826 ms (6492402 allocations: 280.30 MiB)
# almost_periodic     plane         207.373 ms (11530615 allocations: 469.35 MiB)
# almost_periodic     projective    1.091 s (11530616 allocations: 567.03 MiB)
# near_attractor      plane         1.411 s (95942828 allocations: 6.25 GiB)
# near_attractor      projective    1.529 s (96442379 allocations: 6.30 GiB)
# Almost no changes, so we will remove the output type annotations.

# After changing distance for is_nearby
# Criterion           Metric        Benchmark
# escape_time         plane         167.518 ms (6492401 allocations: 241.22 MiB)
# escape_time         projective    665.109 ms (6591444 allocations: 284.09 MiB)
# almost_periodic     plane         207.751 ms (11530615 allocations: 469.35 MiB)
# almost_periodic     projective    1.083 s (11530616 allocations: 567.03 MiB)
# near_attractor      plane         358.019 ms (7050615 allocations: 283.80 MiB)
# near_attractor      projective    573.090 ms (7050616 allocations: 322.89 MiB)
# Why do we get significant speed increases and less allocations?

# After making Attractor mutable
# Criterion           Metric        Benchmark
# escape_time         plane         152.065 ms (6131575 allocations: 202.67 MiB)
# escape_time         projective    663.850 ms (6131576 allocations: 241.76 MiB)
# almost_periodic     plane         209.421 ms (11530682 allocations: 469.35 MiB)
# almost_periodic     projective    1.094 s (11530683 allocations: 567.04 MiB)
# near_attractor      plane         331.953 ms (6410682 allocations: 215.45 MiB)
# near_attractor      projective    547.363 ms (6410683 allocations: 254.54 MiB)
# Makes it faster and uses less allocations