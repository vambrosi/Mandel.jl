module Refactor

export Viewer, Attractor, find_attractors, critical_points

using GLMakie, Symbolics, StaticArraysCore, LinearAlgebra, Polynomials
using GLMakie.Colors
using GLMakie.Colors: LCHab
import Dates

GLMakie.activate!(title = "MandelMakie")

# --------------------------------------------------------------------------------------- #
# Complex Plane Definitions
# --------------------------------------------------------------------------------------- #

const ∞ = Inf + Inf * im

function distance(z::Number, w::Number)
    !isfinite(w) && return 1 / abs(z)
    return abs(z - w)
end

function orbit(f::Function, z::Number, c::Number, iterates::Integer)
    zs = Vector{ComplexF64}(undef, iterates + 1)
    zs[1] = z

    for i = 2:(iterates+1)
        zs[i] = f(zs[i-1], c)
    end

    return zs
end

# --------------------------------------------------------------------------------------- #
# Complex Projective Line Definitions
# --------------------------------------------------------------------------------------- #

const Point = SVector{2,ComplexF64}

function Base.show(io::IO, pt::Point)
    z = pt[1] / pt[2]

    isfinite(z) ? print(io, z) : print(io, "∞")
end

Base.show(io::IO, ::MIME"text/plain", pt::Point) =
    print(io, "Point on the complex projective line:\n   ", pt)

Base.convert(::Type{Point}, z::ComplexF64) =
    isfinite(z) ? normalize(Point(z, 1)) : Point(1, 0)

function distance(pt1::Point, pt2::Point)
    return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

function extend_function(f::Function)
    @variables u, v
    frac = simplify(f(u / v))
    value = Symbolics.value(frac)

    if value isa Number
        h = let
            pt = convert(Point, value)
            _ -> pt
        end
    else
        num, den = Symbolics.arguments(value)
        fu = build_function(num, u, v, expression = Val{false})
        fv = build_function(den, u, v, expression = Val{false})

        h = z -> normalize(Point(fu(z...), fv(z...)))
    end

    g(z::ComplexF64) = f(z)
    g(z::Point) = h(z)
    return g
end

function extend_family(f::Function)
    @variables u, v, a, b
    frac = simplify(f(u / v, a / b))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    fu = build_function(num, u, v, a, b, expression = Val{false})
    fv = build_function(den, u, v, a, b, expression = Val{false})

    g(z::Number, c::Number) = f(z, c)
    g(z::Point, c::Number) = normalize(Point(fu(z..., c, 1), fv(z..., c, 1)))
    g(z::Point, c::Point) = normalize(Point(fu(z..., c...), fv(z..., c...)))

    return g
end

# --------------------------------------------------------------------------------------- #
# Dynamics
# --------------------------------------------------------------------------------------- #

struct DynamicalSystem
    map::Function
    critical_point::Function

    function DynamicalSystem(f::Function, critical_point::Function)
        if hasmethod(f, ComplexF64) && !hasmethod(f, Tuple{ComplexF64,ComplexF64})
            h = (z, c) -> f(z)
        else
            h = (z, c) -> f(z, c)
        end

        return new(extend_family(h), extend_function(critical_point))
    end
end

function DynamicalSystem(f::Function, c::Number)
    c = convert(ComplexF64, c)
    g = let c = c
        _ -> c
    end

    return DynamicalSystem(f, g)
end

const PointLike = Union{ComplexF64,Point}
Base.convert(::Type{PointLike}, z::Number) = ComplexF64(z)
Base.convert(::Type{ComplexF64}, z::Point) = to_complex_plane(z)

struct Attractor{T}
    cycle::Union{T,Vector{T}}
    period::Int
    multiplier::ComplexF64
    power::Int
    palette::Vector{RGB{Float64}}
end

const twilight_RGB = convert.(RGB{Float64}, cgrad(:twilight))
Attractor(c::T, m::Number, p::Integer) where {T} = Attractor{T}(c, 1, m, p, twilight_RGB)
Attractor(c::Vector{T}, m::Number, p::Integer) where {T} =
    Attractor{T}(c, length(c), m, p, twilight_RGB)

function Attractor(
    c::Vector{T},
    m::Number,
    p::Integer,
    index_color::Int,
    n_colors::Int,
) where {T}
    palette = create_gradient(index_color, n_colors)
    Attractor{T}(c, length(c), m, p, palette)
end

function Base.show(io::IO, attractor::Attractor{T}) where {T<:PointLike}
    display_string = "... ↦ "

    for z in attractor.cycle
        display_string *= string(z) * " ↦ "
    end

    display_string *= "..."

    print(io, display_string)
end

Base.show(io::IO, ::MIME"text/plain", attractor::Attractor{ComplexF64}) =
    print(io, "Attracting cycle in the complex plane:\n   ", attractor)

Base.show(io::IO, ::MIME"text/plain", attractor::Attractor{Point}) =
    print(io, "Attracting cycle in the complex projective line:\n   ", attractor)

function Base.convert(::Type{Attractor{ComplexF64}}, attractor::Attractor{Point})
    cycle = convert(Vector{ComplexF64}, attractor.cycle)
    return Attractor(
        cycle,
        attractor.period,
        attractor.multiplier,
        attractor.power,
        attractor.palette,
    )
end

const MaybeAttractor = Union{Nothing,Attractor{ComplexF64},Attractor{Point}}

const empty_attractor = nothing

distance(z, attractor::Attractor{T}) where {T<:PointLike} = distance(z, attractor.cycle)
distance(cycle1::Vector{T}, cycle2::Vector{T}) where {T<:PointLike} =
    minimum([distance(z, w) for z in cycle1, w in cycle2])
distance(z::T, cycle::Vector{T}) where {T<:PointLike} =
    minimum([distance(z, w) for w in cycle])

function convergence_time(
    f::Function,
    z::T,
    c::T,
    attractors::Vector{Attractor{T}},
    ε::Float64,
    max_iterations::Int,
)::Tuple{Int,Float64,MaybeAttractor} where {T<:PointLike}
    for iteration = 0:max_iterations
        for attractor in attractors
            d = distance(z, attractor)
            if d < ε
                return iteration, d, attractor
            end
        end

        z = f(z, c)
    end

    return max_iterations + 1, ε, empty_attractor
end

struct CloseBy
    iterations::Int
    z::PointLike
    w::PointLike
end

const MaybeCloseBy = Union{Nothing,CloseBy}

function period_multiple_apart(f, z, c, ε, max_iterations)::MaybeCloseBy
    slow = fast = z

    for n = 1:max_iterations
        slow = f(slow, c)
        fast = f(f(fast, c), c)

        distance(slow, fast) <= ε && return CloseBy(n, slow, fast)
    end

    return nothing
end

function iterate_until_close(f, z, w, c, ε, max_iterations)::MaybeCloseBy
    for n = 1:max_iterations
        z = f(z, c)

        distance(z, w) < ε && return CloseBy(n, z, w)
    end

    return nothing
end

function iterate_both_until_close(f, z, w, c, ε, max_iterations)::MaybeCloseBy
    for n = 1:max_iterations
        z = f(z, c)
        w = f(w, c)

        distance(z, w) < ε && return CloseBy(n, z, w)
    end

    return nothing
end

struct OrbitData
    preperiod::Int
    period::Int
end

const MaybeOrbitData = Union{Nothing,OrbitData}

function multiplier(f, z0, c, _, ε, max_iterations)::MaybeOrbitData
    points = period_multiple_apart(f, z0, c, ε, max_iterations)
    isnothing(points) && return nothing

    points = iterate_until_close(f, points.z, points.w, c, ε, max_iterations)
    isnothing(points) && return nothing

    period = points.iterations
    points = iterate_until_close(f, z0, points.w, c, ε, max_iterations)
    isnothing(points) && return nothing

    preperiod = points.iterations
    return OrbitData(preperiod, period)
end

# --------------------------------------------------------------------------------------- #
# Finding Attractors
# --------------------------------------------------------------------------------------- #

function coeffs(polynomial, z)
    polynomial = expand(polynomial)
    d = Symbolics.degree(polynomial)

    coefficients = [Symbolics.coeff(polynomial, z^d) for d = 1:d]
    prepend!(coefficients, Symbolics.substitute(polynomial, Dict(z => 0)))
    return convert.(ComplexF64, coefficients)
end

function critical_points(func, parameter)
    @variables z, u, v, c
    f = func(z, c)

    # Hack to write rational function as a ratio of polynomials
    # (Not sure it works on all cases.)
    f =
        f |>
        (x -> substitute(x, Dict(z => u / v))) |>
        simplify |>
        (x -> substitute(x, Dict(u => z, v => 1)))

    # Get numerator and denominator polynomials
    p, q = f |> Symbolics.value |> Symbolics.arguments

    parameter = convert(ComplexF64, parameter)
    p = Polynomial(coeffs(substitute(p, Dict(c => parameter)), z))
    q = Polynomial(coeffs(substitute(q, Dict(c => parameter)), z))

    num = Polynomials.derivative(p) * q - p * Polynomials.derivative(q)
    points = unique!(convert.(ComplexF64, roots(num)))

    d = max(Polynomials.degree(p), Polynomials.degree(q))
    length(points) < 2 * d - 2 && pushfirst!(points, ∞)
    return points
end

function attracting_cycle(f, z::T, c, ε, max_iterations) where {T<:PointLike}
    points = period_multiple_apart(f, z, c, ε, max_iterations)
    isnothing(points) && return T[]

    orbiter = points.z
    reference = points.w

    orbit = [orbiter]
    for _ = 1:max_iterations
        orbiter = f(orbiter, c)
        distance(orbiter, reference) <= ε && return orbit
        push!(orbit, orbiter)
    end

    return orbit
end

function find_attractors(f::Function; projective::Bool = false)
    hasmethod(f, ComplexF64) || throw("If it is a family of functions, input a parameter.")

    h = extend_family((z, c) -> f(z))

    return find_attractors(h, 0.0im; projective = projective)
end

function find_attractors(f::Function, c::Number; projective::Bool = false)
    # If it is not a family ignores the parameter
    hasmethod(f, Tuple{ComplexF64,ComplexF64}) || return find_attractors(f)

    c = convert(ComplexF64, c)

    # Extend the function, in case it was not extended previously
    g = extend_family(f)

    # Iterate critical points to find all attracting cycles
    crit_pts = critical_points(g, c)

    orbits = Vector{Point}[]
    for z in crit_pts
        pt = convert(Point, z)
        orbit = attracting_cycle(g, pt, c, 1e-4, 1000)

        if length(orbit) > 0
            is_repeat = any([distance(orbit, o) < 2e-4 for o in orbits])
            !is_repeat && push!(orbits, orbit)
        end
    end

    projective || (orbits = convert.(Vector{ComplexF64}, orbits))

    n = length(orbits)
    n == 1 && return [Attractor(orbits[1], 0.0im, 0)]
    return [Attractor(orbit, 0.0im, 0, i, n) for (i, orbit) in enumerate(orbits)]
end

# --------------------------------------------------------------------------------------- #
# Coloring Algorithms
# --------------------------------------------------------------------------------------- #

const DEFAULT_COLOR_LEVEL = 0.5
const DEFAULT_COLOR = twilight_RGB[255]

function to_color(approach::OrbitData)
    depth = mod(approach.preperiod / approach.period / 64.0, 1.0)
    return twilight_RGB[round(Int, 509 * depth + 1)]
end

to_color(::Integer, ::Number, ::Nothing) = DEFAULT_COLOR

function to_color(iterations::Integer, d::Number, attractor::Attractor)
    depth = mod(
        (iterations / attractor.period + 1.0 - log(attractor.power, -log(abs(d)))) / 64.0,
        1.0,
    )
    return attractor.palette[round(Int, 509 * depth + 1)]
end

to_color(::Nothing) = DEFAULT_COLOR

escape_time(f, z, c, a, ε, N) = to_color(convergence_time(f, z, c, a, ε, N)...)
convergence_color(f, z, c, a, ε, N) = to_color(multiplier(f, z, c, a, ε, N))

struct ColoringData{T}
    method::Function
    attractors::Vector{Attractor{T}}
    update_attractors::Bool
end

attractor_type(::ColoringData{T}) where {T} = T
attractor_type(is_projective::Bool) = is_projective ? Point : ComplexF64

const base_hue_chroma = [
    #RGB CA736C 202 115 108 red
    (38.35791098655309, 29.49442824523494),
    #RGB 5794D0 87 148 208 blue
    (36.822489073365745, 266.89189631383687),
    #RGB 47A477 71 164 119 green
    (41.31205868565174, 158.3817555910295),
    #RGB 8D9741 141 151 65 yellow
    (46.2028504263141, 110.51854840517339),
    #RGB 00A2AF 0 162 175 cyan
    (34.439876883846594, 209.09162086774987),
    #RGB BC73A4 188 115 164 magenta
    (38.49553456026353, 338.47504267487597),
    #RGB BA823A 186 130 58 orange
    (48.58788683488929, 72.5301368954557),
    #RGB 9481CC 148 129 204 purple
    (43.17999997459473, 302.8931595387337),
]

function create_gradient(color_index, n_colors)
    if n_colors == 1
        return twilight_RGB
    elseif n_colors < 9
        chroma, hue = base_hue_chroma[color_index]
    else
        chroma, hue = 40, mod(360 * (color_index - 1) / n + 30, 360)
    end

    gradient = [LCHab(20 * cospi(t) + 50, chroma, hue) for t in range(0.0, 2.0, 510)]
    return convert.(RGB{Float64}, gradient)
end

# --------------------------------------------------------------------------------------- #
# Views of Mandelbrot and Julia Sets
# --------------------------------------------------------------------------------------- #

abstract type View end

mutable struct Options
    convergence_radius::Float64
    max_iterations::Int
    orbit_length::Int
    critical_length::Int
    compact_view::Bool
    is_family::Bool
end

mutable struct MandelView <: View
    center::ComplexF64
    diameter::Float64
    pixels::Int

    init_center::ComplexF64
    init_diameter::Float64

    color_levels::Observable{Matrix{RGB{Float64}}}
    points::Observable{Vector{ComplexF64}}
    marks::Observable{Vector{ComplexF64}}

    coloring_data::ColoringData

    function MandelView(; center, diameter, pixels, coloring_data)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        color_levels = zeros(Float64, pixels, pixels)
        points = ComplexF64[center]
        marks = ComplexF64[]

        return new(
            center,
            diameter,
            pixels,
            center,
            diameter,
            color_levels,
            points,
            marks,
            coloring_data,
        )
    end
end

mutable struct JuliaView <: View
    center::ComplexF64
    diameter::Float64
    parameter::ComplexF64
    pixels::Int

    init_center::ComplexF64
    init_diameter::Float64

    color_levels::Observable{Matrix{RGB{Float64}}}
    points::Observable{Vector{ComplexF64}}
    marks::Observable{Vector{ComplexF64}}

    coloring_data::ColoringData

    function JuliaView(; center, diameter, parameter, pixels, coloring_data)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        color_levels = zeros(Float64, pixels, pixels)
        points = ComplexF64[center]
        marks = ComplexF64[]

        return new(
            center,
            diameter,
            parameter,
            pixels,
            center,
            diameter,
            color_levels,
            points,
            marks,
            coloring_data,
        )
    end
end

# --------------------------------------------------------------------------------------- #
# Change of Coordinates between Complex Plane, Projective Line, and Pixel Space
# --------------------------------------------------------------------------------------- #

function to_complex_plane(view::View, pixel_vector)
    x, y = pixel_vector[1], pixel_vector[2]

    upp = view.diameter / view.pixels
    a = (x + 0.5 - view.pixels / 2) * upp
    b = (y + 0.5 - view.pixels / 2) * upp

    return view.center + complex(a, b)
end

function to_complex_plane(pt::Point)
    z = pt[1] / pt[2]

    !isfinite(z) && return ∞
    return z
end

function to_pixel_space(view::View, z::Number)
    ppu = view.pixels / view.diameter
    w = z - view.center

    x = real(w) * ppu - 0.5 + view.pixels / 2
    y = imag(w) * ppu - 0.5 + view.pixels / 2

    return [x, y]
end

function to_pixel_space(view::View, zs::Vector{<:Number})
    ppu = view.pixels / view.diameter
    xs = Float64[]
    ys = Float64[]

    for z in zs
        w = z - view.center
        push!(xs, real(w) * ppu - 0.5 + view.pixels / 2)
        push!(ys, imag(w) * ppu - 0.5 + view.pixels / 2)
    end

    return xs, ys
end

# --------------------------------------------------------------------------------------- #
# Updating Plots
# --------------------------------------------------------------------------------------- #

function pick_orbit!(
    julia::JuliaView,
    d_system::DynamicalSystem,
    options::Options,
    z::Number,
)
    julia.points[] = orbit(d_system.map, z, julia.parameter, options.orbit_length - 1)
    return julia
end

function corner_and_step(view::View)
    step = view.diameter / view.pixels

    corner_real = real(view.center) - view.diameter / 2 + 0.5 * step
    corner_imag = imag(view.center) - view.diameter / 2 + 0.5 * step
    corner = complex(corner_real, corner_imag)

    return corner, step
end

function mandel_slice!(array, j, f, crit, corner, step, pxs, coloring_data, options)
    @inbounds for i = 1:pxs
        c = convert(attractor_type(coloring_data), corner + step * complex(i, j))

        array[i, j] = coloring_data.method(
            f,
            crit(c),
            c,
            coloring_data.attractors,
            options.convergence_radius,
            options.max_iterations,
        )
    end
    return nothing
end

function update_grid!(
    view::MandelView,
    d_system::DynamicalSystem,
    corner::ComplexF64,
    step::Float64,
    options::Options,
)
    futures = Vector{Task}(undef, view.pixels)

    @inbounds for j = 1:view.pixels
        futures[j] = Threads.@spawn mandel_slice!(
            view.color_levels[],
            j,
            d_system.map,
            d_system.critical_point,
            corner,
            step,
            view.pixels,
            view.coloring_data,
            options,
        )
    end

    wait.(futures)
    return nothing
end

function julia_slice!(array, j, f, c, corner, step, pxs, coloring_data, options)
    @inbounds for i = 1:pxs
        z = convert(attractor_type(coloring_data), corner + step * complex(i, j))

        array[i, j] = coloring_data.method(
            f,
            z,
            c,
            coloring_data.attractors,
            options.convergence_radius,
            options.max_iterations,
        )
    end
    return nothing
end

function update_grid!(
    view::JuliaView,
    d_system::DynamicalSystem,
    corner::ComplexF64,
    step::Float64,
    options::Options,
)
    parameter = convert(attractor_type(view.coloring_data), view.parameter)
    futures = Vector{Task}(undef, view.pixels)

    @inbounds for j = 1:view.pixels
        futures[j] = Threads.@spawn julia_slice!(
            view.color_levels[],
            j,
            d_system.map,
            parameter,
            corner,
            step,
            view.pixels,
            view.coloring_data,
            options,
        )
    end

    wait.(futures)
    return nothing
end

function update_view!(view::View, d_system::DynamicalSystem, options::Options)
    corner, step = corner_and_step(view)
    update_grid!(view, d_system, corner, step, options)
    notify(view.color_levels)
    notify(view.points)
    notify(view.marks)

    return view
end

function pick_parameter!(
    julia::JuliaView,
    mandel::MandelView,
    d_system::DynamicalSystem,
    options::Options,
    point,
)
    julia.parameter = point

    julia.marks[] = orbit(
        d_system.map,
        d_system.critical_point(julia.parameter),
        julia.parameter,
        options.critical_length - 1,
    )

    mandel.points[] = [julia.parameter]

    coloring_data = julia.coloring_data
    if coloring_data.update_attractors
        T = attractor_type(coloring_data)
        attractor_list =
            find_attractors(d_system.map, julia.parameter, projective = (T == Point))
        julia.coloring_data = ColoringData{T}(
            coloring_data.method,
            attractor_list,
            coloring_data.update_attractors,
        )
    end

    update_view!(julia, d_system, options)
    pick_orbit!(julia, d_system, options, julia.points[][begin])
    return julia
end

# --------------------------------------------------------------------------------------- #
# User Interface
# --------------------------------------------------------------------------------------- #

struct Frame
    axis::Axis
    view::Ref{Any}
    events::Dict{Symbol,Any}
end

function translate_axis!(axis, z)
    translate!(axis.scene, 0, 0, z)
    translate!(axis.elements[:background], 0, 0, z - 1)
    translate!(axis.elements[:xoppositeline], 0, 0, z + 1)
    translate!(axis.elements[:yoppositeline], 0, 0, z + 1)
    translate!(axis.xaxis.elements[:axisline], 0, 0, z + 1)
    translate!(axis.yaxis.elements[:axisline], 0, 0, z + 1)
end

function create_frames!(figure, options, mandel::Nothing, julia)
    julia_limits = (0.5, 0.5 + julia.pixels, 0.5, 0.5 + julia.pixels)
    axis = Axis(figure[1, 1][1, 1], aspect = 1, limits = julia_limits)
    translate_axis!(axis, 0)

    hidedecorations!(axis)
    deregister_interaction!(axis, :rectanglezoom)
    deregister_interaction!(axis, :limitreset)
    deregister_interaction!(axis, :dragpan)

    return nothing, Frame(axis, julia, Dict())
end

function create_frames!(figure, options, mandel::MandelView, julia)
    mandel_limits = (0.5, 0.5 + mandel.pixels, 0.5, 0.5 + mandel.pixels)
    julia_limits = (0.5, 0.5 + julia.pixels, 0.5, 0.5 + julia.pixels)

    left_axis = Axis(figure[1, 1][1, 1], aspect = 1, limits = mandel_limits)
    right_axis = Axis(figure[1, 1][1, 1], aspect = 1, limits = julia_limits)
    translate_axis!(left_axis, 10)
    translate_axis!(right_axis, 0)

    if options.compact_view
        left_axis.width = Relative(0.3)
        left_axis.height = Relative(0.3)
        left_axis.valign = 0.03
        left_axis.halign = 0.03
    end

    for axis in [left_axis, right_axis]
        hidedecorations!(axis)
        deregister_interaction!(axis, :rectanglezoom)
        deregister_interaction!(axis, :limitreset)
        deregister_interaction!(axis, :dragpan)
    end

    return Frame(left_axis, mandel, Dict()), Frame(right_axis, julia, Dict())
end

function create_plot!(frame::Frame)
    empty!(frame.axis)
    view = frame.view[]

    heatmap!(
        frame.axis,
        view.color_levels,
        colormap = :twilight,
        colorrange = (0.0, 1.0),
        inspectable = false,
        inspector_label = (self, i, p) -> let
            z = to_complex_plane(view, p)
            "x: $(real(z))\ny: $(imag(z))"
        end,
    )

    point_vectors = lift(view.points) do zs
        xs, ys = to_pixel_space(view, zs)
        return Point2f.(xs, ys)
    end

    lines!(frame.axis, point_vectors, color = (:red, 0.5), inspectable = false)

    scatter!(
        frame.axis,
        point_vectors,
        color = (:red, 1.0),
        inspector_label = (self, i, p) -> let
            z = to_complex_plane(view, p)
            "x: $(real(z))\ny: $(imag(z))"
        end,
    )

    mark_vectors = lift(view.marks) do zs
        xs, ys = to_pixel_space(view, zs)
        return Point2f.(xs, ys)
    end

    lines!(frame.axis, mark_vectors, color = (:blue, 0.5), inspectable = false)

    scatter!(
        frame.axis,
        mark_vectors,
        color = (:blue, 1.0),
        inspector_label = (self, i, p) -> let
            z = to_complex_plane(view, p)
            "x: $(real(z))\ny: $(imag(z))"
        end,
    )

    return frame
end

function zoom!(frame::Frame, d_system::DynamicalSystem, options::Options)
    frame.events[:is_zooming] = true
    view = frame.view[]

    x_min, x_max, _, _ = frame.axis.limits[]
    original_size = x_max - x_min
    x1_min, x1_max = frame.axis.xaxis.attributes.limits[]
    current_size = x1_max - x1_min

    scale = current_size / original_size
    vector = mouseposition(frame.axis.scene)
    z = to_complex_plane(view, vector)

    view.diameter *= scale
    view.center = scale * view.center + (1 - scale) * z
    update_view!(view, d_system, options)
    reset_limits!(frame.axis)
    frame.events[:is_zooming] = false

    return frame
end

function save_view(filename::String, view::View)
    pxs = view.pixels
    fig = Figure(figure_padding = 0, size = (pxs, pxs))
    ax = Axis(fig[1, 1], aspect = AxisAspect(1))

    hidedecorations!(ax)
    hidespines!(ax)

    heatmap!(ax, view.color_levels[], colormap = :twilight, colorrange = (0.0, 1.0))

    Makie.save(filename, fig)
end

function add_frame_events!(
    frame::Frame,
    topframe::Frame,
    d_system::DynamicalSystem,
    options::Options,
    julia::JuliaView,
)
    axis = frame.axis
    scene = axis.scene

    # Mouse Events
    dragging = false
    dragstart = Point2f(0.0)
    dragend = Point2f(0.0)

    to_world_at_start = z -> to_world(scene, z)

    is_zooming = false
    zooming = Timer(_ -> zoom!(frame, d_system, options), 0.1)

    is_topframe = frame == topframe
    z_level = is_topframe ? 10 : 0

    on(events(scene).mousebutton) do event
        is_zooming && return Consume(false)

        mp = mouseposition_px(scene)
        view = frame.view[]

        if event.button == Mouse.right
            if event.action == Mouse.press &&
               is_mouseinside(axis) &&
               (is_topframe || !is_mouseinside(topframe.axis))
                dragging = true
                to_world_at_start = z -> to_world(scene, z)
                dragstart = to_world_at_start(mp)

            elseif event.action == Mouse.release && dragging
                dragend = to_world_at_start(mp)
                view.center +=
                    to_complex_plane(view, dragstart) - to_complex_plane(view, dragend)
                translate!(scene, 0, 0, z_level)
                update_view!(view, d_system, options)
                reset_limits!(axis)
                dragging = false
            end
        end

        if event.button == Mouse.left
            if event.action == Mouse.press &&
               is_mouseinside(axis) &&
               (is_topframe || !is_mouseinside(topframe.axis))

                if ispressed(scene, Keyboard.left_control | Keyboard.right_control)
                    view.center = view.init_center
                    view.diameter = view.init_diameter

                    translate!(scene, 0, 0, z_level)
                    update_view!(view, d_system, options)
                    reset_limits!(axis)

                    return Consume(true)
                end

                point = to_complex_plane(view, to_world_at_start(mp))
                if view isa MandelView
                    pick_parameter!(julia, view, d_system, options, point)
                elseif view isa JuliaView
                    pick_orbit!(view, d_system, options, point)
                end
            end
        end

        return Consume(false)
    end

    on(events(scene).mouseposition) do event
        if dragging
            mp = mouseposition(scene)
            translate!(scene, mp - dragstart..., z_level)
        end
    end

    on(events(scene).scroll, priority = 100) do event
        if is_mouseinside(axis) &&
           !is_zooming &&
           (is_topframe || !is_mouseinside(topframe.axis))
            close(zooming)
            zooming = Timer(_ -> zoom!(frame, d_system, options), 0.1)
        end
    end

    frame.events[:dragging] = dragging
    frame.events[:dragstart] = dragstart
    frame.events[:dragend] = dragstart
    frame.events[:is_zooming] = is_zooming
    frame.events[:zooming] = zooming

    # Keyboard Events
    on(events(scene).keyboardbutton) do event
        is_zooming && return Consume(false)
        view = frame.view[]

        if ispressed(scene, (Keyboard.left_control | Keyboard.right_control) & Keyboard.s)
            if is_mouseinside(axis) &&
               !is_zooming &&
               (is_topframe || !is_mouseinside(topframe.axis))

                format = Dates.dateformat"yyyy-mm-ddTHH.MM.SS"
                time = string(Dates.format(Dates.now(), format))
                !isdir("imgs") && mkdir("imgs")
                save_view("imgs/" * time * ".png", view)
                return Consume(true)
            end
        end
    end
end

function add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)
    layout = GridLayout(figure[2, 1], tellwidth = false)
    button_shift = options.is_family ? 2 : 0

    labels = Dict(
        :max_iter => Label(layout[1, button_shift+1], "Maximum\nIterations:"),
        :orbit_len => Label(layout[1, button_shift+3], "Orbit\nLength:"),
        :critical_length =>
            Label(layout[1, button_shift+5], "Critical Point\nOrbit Length:"),
        :convergence_radius => Label(layout[1, button_shift+7], "Convergence\nRadius:"),
    )

    inputs = Dict{Symbol,Any}(
        :max_iter => Textbox(
            layout[1, button_shift+2],
            width = 60,
            placeholder = string(options.max_iterations),
            validator = Int,
        ),
        :orbit_len => Textbox(
            layout[1, button_shift+4],
            width = 60,
            placeholder = string(options.orbit_length),
            validator = Int,
        ),
        :critical_length => Textbox(
            layout[1, button_shift+6],
            width = 60,
            placeholder = string(options.orbit_length),
            validator = Int,
        ),
        :convergence_radius => Textbox(
            layout[1, button_shift+8],
            width = 60,
            placeholder = string(options.convergence_radius),
            validator = Float64,
        ),
    )

    on(inputs[:max_iter].stored_string) do s
        options.max_iterations = parse(Int, s)
        update_view!(julia, d_system, options)
        options.is_family && update_view!(mandel, d_system, options)
    end

    on(inputs[:orbit_len].stored_string) do s
        options.orbit_length = parse(Int, s)
        pick_orbit!(julia, d_system, options, julia.points[][begin])
    end

    on(inputs[:critical_length].stored_string) do s
        options.critical_length = parse(Int, s)
        julia.marks[] = orbit(
            d_system.map,
            d_system.critical_point(julia.parameter),
            julia.parameter,
            options.critical_length - 1,
        )
    end

    on(inputs[:convergence_radius].stored_string) do s
        options.convergence_radius = parse(Float64, s)
        update_view!(julia, d_system, options)
        options.is_family && update_view!(mandel, d_system, options)
    end

    if options.is_family
        inputs[:switch_layout] = Button(layout[1, 1], label = "↰", halign = :left)
        inputs[:switch_positions] = Button(layout[1, 2], label = "↔", halign = :left)

        on(inputs[:switch_layout].clicks, priority = 300) do event
            if options.compact_view
                left_frame.axis.width = nothing
                left_frame.axis.height = nothing
                left_frame.axis.valign = :center
                left_frame.axis.halign = :center

                figure[1, 1][1, 1] = left_frame.axis
                figure[1, 1][1, 2] = right_frame.axis

                inputs[:switch_layout].label = "↳"
            else
                left_frame.axis.width = Relative(0.3)
                left_frame.axis.height = Relative(0.3)
                left_frame.axis.valign = 0.01
                left_frame.axis.halign = 0.01

                figure[1, 1][1, 1] = left_frame.axis
                figure[1, 1][1, 1] = right_frame.axis

                trim!(contents(figure[1, 1])...)
                inputs[:switch_layout].label = "↰"
            end

            options.compact_view = !options.compact_view
        end

        on(inputs[:switch_positions].clicks, priority = 300) do event
            left_frame.view[], right_frame.view[] = right_frame.view[], left_frame.view[]
            create_plot!(left_frame)
            create_plot!(right_frame)
        end
    end

    return inputs
end

make_tuple(x) = (x, x)
make_tuple(x::Tuple) = x

function fix_criteria((mandel, julia))
    # Test if it is one of the valid options
    options = [:escape_time, :near_attractor, :almost_periodic]
    mandel in options || throw("Invalid Mandelbrot set coloring method")
    julia in options || throw("Invalid Julia set coloring method")

    # Mandelbrot can't use attractors so change to default
    mandel == :near_attractor && (mandel = :almost_periodic)
    return mandel, julia
end

function get_coloring_data(map, c, convergence_criterion, projective_metric)
    if convergence_criterion == :near_attractor
        method = escape_time
        attractors = find_attractors(map, c, projective = projective_metric)
        update_attractors = true
    elseif convergence_criterion == :escape_time
        method = escape_time
        attractors = [Attractor(∞, 0, 2)]
        update_attractors = false
    elseif convergence_criterion == :almost_periodic
        method = convergence_color
        attractors = Attractor{attractor_type(projective_metric)}[]
        update_attractors = false
    else
        throw("Invalid `convergence_criterion`")
    end

    return ColoringData(method, attractors, update_attractors)
end

"""
	Viewer(f; <keyword arguments>)

Create a `Viewer` that plots the Mandelbrot and the Julia sets associated with the \
function `f`. `f` can have one or two inputs, the second being the parameter. \
By default, the `Viewer` opens as a separate window.

# Examples
```julia-repl
julia> Viewer((z, c) -> z^2 + c, mandel_center=-0.5)
```
```julia-repl
julia> f(z, λ) = z^2 + λ / z^2
julia> crit(λ) = λ^(1/4)
julia> Viewer(f; crit=crit, mandel_diameter=1.0)
```

# Shortcuts

- `Right Click Drag`: Pans view.
- `Mouse Scroll`: Zooms in or out.
- `Left Click` (parameter space): Chooses parameter.
- `Left Click` (dynamical space): Chooses orbit initial point.
- `Ctrl + Left Click`: Resets view to initial `center` and `diameter`.
- `Ctrl + S`: Saves view that is under the mouse pointer (files in `./imgs`).

# Arguments
- `crit = 0.0im`: Function that gives a critical point for each parameter. Used to \
to plot the Mandelbrot set. If it is a constant function, you can just input the constant \
directly.

- `c = 0.0im`: Initial parameter used to plot the Julia set.

- `mandel_center = 0.0im`: Initial center of the Mandelbrot plot.

- `mandel_diameter = 4.0`: Initial diameter of the Mandelbrot plot.

- `julia_center = 0.0im`: Initial center of the Julia plot.

- `julia_diameter = 4.0`: Initial diameter of the Julia plot.

- `grid_width = 800`: Width (and height) of the grid of complex numbers used to plot sets.

- `compact_view = true`: If 'true' one of the plots is show as an inset plot, if `false` \
they are shown side-by-side.

# Coloring Method Options
For the options below, if you want to set different values for the Mandelbrot and Julia \
set views, set the option to be a tuple with the respective values.

- `convergence_criterion = :escape_time`: Chooses the coloring method for both plots. \
The options are `:escape_time`, `:near_attractor`, `:almost_periodic`. More details below.

- `projective_metric = false`: Determines which metric will be used to determine \
distances, the complex plane metric, or the metric on the projective line. The distance \
between ∞ and a finite point in the plane metric is the inverse of its absolute value.

# Covergence Criteria

- `:escape_time` (default): computes how fast a point approaches ∞ (in the plane metric);
- `:almost_periodic`: finds the attracting cycle for each point separately \
using Floyd's cycle-finding algorithm;
- `:near_attractor`: computes all attracting cycles in advance, and then computes how \
fast each point converges to one of those attractors (uses different color gradients \
for each attractor). This option is only available for the Julia set, and it will \
default to `:almost_periodic` in the Mandelbrot set case.

"""
struct Viewer
    d_system::DynamicalSystem
    options::Options

    figure::Figure
    left_frame::Union{Nothing,Frame}
    right_frame::Frame

    mandel::Union{Nothing,MandelView}
    julia::JuliaView
    inputs::Dict{Symbol,Any}

    function Viewer(
        f;
        crit = 0.0im,
        c = 0.0im,
        mandel_center = 0.0im,
        mandel_diameter = 4.0,
        julia_center = 0.0im,
        julia_diameter = 4.0,
        compact_view = true,
        grid_width = 800,
        convergence_criterion = :escape_time,
        projective_metric = false,
    )
        # Put options in standard form
        projective_metrics = make_tuple(projective_metric)
        convergence_criteria = fix_criteria(make_tuple(convergence_criterion))

        # Check if it is a family of maps
        is_family = hasmethod(f, Tuple{ComplexF64,ComplexF64})

        # Create Viewer Data
        d_system = DynamicalSystem(f, crit)
        options = Options(1e-3, 200, 1, 1, compact_view, is_family)
        figure = Figure(size = (750, 650))

        mandel_coloring = get_coloring_data(
            d_system.map,
            c,
            convergence_criteria[1],
            projective_metrics[1],
        )
        julia_coloring = get_coloring_data(
            d_system.map,
            c,
            convergence_criteria[2],
            projective_metrics[2],
        )

        mandel =
            !is_family ? nothing :
            MandelView(
                center = mandel_center,
                diameter = mandel_diameter,
                pixels = grid_width,
                coloring_data = mandel_coloring,
            )

        julia = JuliaView(
            center = julia_center,
            diameter = julia_diameter,
            parameter = c,
            pixels = grid_width,
            coloring_data = julia_coloring,
        )

        julia.marks[] = [d_system.critical_point(julia.parameter)]
        pick_orbit!(julia, d_system, options, julia.points[][begin])

        left_frame, right_frame = create_frames!(figure, options, mandel, julia)

        if is_family
            add_frame_events!(right_frame, left_frame, d_system, options, julia)
            add_frame_events!(left_frame, left_frame, d_system, options, julia)

            create_plot!(left_frame)
            update_view!(mandel, d_system, options)
            mandel.points[] = [julia.parameter]
        else
            add_frame_events!(right_frame, right_frame, d_system, options, julia)
        end

        create_plot!(right_frame)
        update_view!(julia, d_system, options)

        colgap!(content(figure[1, 1]), 10)
        inputs =
            add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)

        return new(
            d_system,
            options,
            figure,
            left_frame,
            right_frame,
            mandel,
            julia,
            inputs,
        )
    end
end

Base.show(io::IO, viewer::Viewer) = display(GLMakie.Screen(), viewer.figure)

end
