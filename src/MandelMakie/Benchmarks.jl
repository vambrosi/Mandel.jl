using BenchmarkTools
using MandelMakie

pick_parameter! = MandelMakie.Refactor.pick_parameter!

function test(c, p)
    f(z, c) = z^2 + c
    viewer = Viewer(
        f,
        mandel_center = -0.5,
        convergence_criterion = c,
        projective_metric = p,
    )

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
    viewer = Viewer(
        f,
        mandel_center = -0.5,
        convergence_criterion = c,
        projective_metric = p,
    )

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
