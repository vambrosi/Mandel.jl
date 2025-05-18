module Refactor

export Viewer, Attractor, get_attractors, get_parameter, critical_points

using GLMakie, Symbolics, StaticArraysCore, LinearAlgebra, Polynomials, HypertextLiteral
include("./Rays.jl")
using GLMakie.Colors
using Crayons

import Dates, Nemo

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

    for i in 2:(iterates+1)
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

    g(z) = f(z)
    g(z::Point) = h(z)
    return g
end

function extend_family(f::Function)
    @variables u, v, a, b
    frac = simplify(f(u / v, a / b))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    f1 = build_function(num, u, v, a, b, expression = Val{false})
    f2 = build_function(den, u, v, a, b, expression = Val{false})

    g(z, c) = f(z, c)
    g(z::Point, c::Number) = normalize(Point(f1(z..., c, 1), f2(z..., c, 1)))
    g(z::Point, c::Point) = normalize(Point(f1(z..., c...), f2(z..., c...)))

    df_symb_00 = expand_derivatives(Differential(u)(simplify(f(u, a))))
    df_symb_0∞ = expand_derivatives(Differential(u)(simplify(1 / f(u, a))))
    df_symb_∞0 = expand_derivatives(Differential(u)(simplify(f(1 / u, a))))
    df_symb_∞∞ = expand_derivatives(Differential(u)(simplify(1 / f(1 / u, a))))

    df_00 = build_function(df_symb_00, u, a, expression = Val{false})
    df_0∞ = build_function(df_symb_0∞, u, a, expression = Val{false})
    df_∞0 = build_function(df_symb_∞0, u, a, expression = Val{false})
    df_∞∞ = build_function(df_symb_∞∞, u, a, expression = Val{false})

    function g(orbit::Vector{ComplexF64}, c::ComplexF64)
        n = length(orbit)
        λ = 1

        for i in 1:n
            z, w = orbit[i], orbit[mod1(i + 1, n)]
            if abs(z) < 1
                λ *= abs(w) < 1 ? df_00(z, c) : df_0∞(z, c)
            else
                λ *= abs(w) < 1 ? df_∞0(z, c) : df_∞∞(z, c)
            end
        end

        return λ
    end

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

        return new(h, critical_point)
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

mutable struct Attractor{T}
    cycle::Union{T,Vector{T}}
    period::Int
    multiplier::Float64
    power::Int
    color::Union{Symbol,LCHab}
    palette::Vector{RGBA{Float64}}
    continuous_coloring::Bool
end

const twilight_RGB = convert.(RGBA{Float64}, cgrad(:twilight))
Attractor(c::T, m::Number, p::Integer) where {T} =
    Attractor{T}(c, 1, m, p, :twilight, twilight_RGB, true)
Attractor(c::Vector{T}, m::Number, p::Integer) where {T} =
    Attractor{T}(c, length(c), m, p, :twilight, twilight_RGB, true)

function Attractor(
    c::Vector{T},
    m::Number,
    p::Integer,
    index_color::Int,
    n_colors::Int,
) where {T}
    color, palette = create_gradient(index_color, n_colors)
    Attractor{T}(c, length(c), m, p, color, palette, true)
end

function Base.show(io::IO, attractor::Attractor{T}; indent::Int = 0) where {T<:PointLike}
    indentation = " "^indent

    if attractor.cycle isa Number
        display_string = indentation * "... ↦ " * string(attractor.cycle) * " ↦ ..."
    else
        display_string = indentation * "... ↦ " * string(attractor.cycle[1]) * " ↦ "

        for z in attractor.cycle[2:end]
            display_string *= "\n" * indentation * "     ↦ " * string(z) * " ↦ "
        end

        display_string *= "..."
    end

    print(io, display_string)
end

function print_base_color(io::IO, color)
    color = convert(RGB24, color).color
    print(io, "base color ")
    print(io, Crayon(foreground = color), "██")
    print(io, Crayon(reset = true))
end

function print_base_color(io::IO, color::Symbol)
    print(io, color)
end

function print_palette(io::IO, palette; size = 15)
    Δ, r = divrem(length(palette) - 1, size)

    position = div(r, 2) + 1
    for i in 0:size
        color = convert(RGB24, palette[position]).color
        print(io, Crayon(foreground = color), "█")
        position += Δ
    end

    print(io, Crayon(reset = true))
end

function Base.show(io::IO, ::MIME"text/plain", attractor::Attractor{T}) where {T<:PointLike}
    space_name = T == ComplexF64 ? "plane" : "projective line"
    print(
        io,
        "Attracting cycle in the complex $space_name: \n  \
        Period     = $(attractor.period) \n  \
        Multiplier = $(attractor.multiplier) \n  \
        Degree     = $(attractor.power) \n  \
        Palette    = $(sprint(print_palette, attractor.palette)) \
        ($(sprint(print_base_color, attractor.color))) \n  \
        Cycle      = " * (attractor.period > 1 ? "\n " : ""),
    )
    show(io, attractor, indent = (attractor.period > 1 ? 3 : 0))
end

function Base.show(io::IO, m::MIME"text/html", attractor::Attractor{T}) where {T<:PointLike}
    space_name = T == ComplexF64 ? "plane" : "projective line"

    if attractor.color isa Symbol
        show(
            io,
            m,
            @htl """<table style="width:550px"> \
                    <tr> \
                        <td colspan=2> Attracting cycle in the complex $space_name: </td> \
                    </tr> \
                    <tr>\
                        <td style="width:20%;"> Period </td>\
                        <td style="width:80%;"> $(attractor.period) </td>\
                    </tr>\
                    <tr>\
                        <td style="width:20%;"> Multiplier </td>\
                        <td style="width:80%;"> $(attractor.multiplier) </td>\
                    </tr>\
                    <tr>\
                        <td style="width:20%;"> Degree </td>\
                        <td style="width:80%;"> $(attractor.power) </td>\
                    </tr>\
                    <tr>
                        <td style="width:20%;"> Color </td>\
                        <td style="width:80%;"> $(attractor.color) </td>\
                    </tr>\
                    <tr><td style="width:20%;"> Cycle  </td> <td style="width:80%;"> """
        )
    else
        color = string(convert(RGB24, attractor.color).color, base = 16)
        show(
            io,
            m,
            @htl """<table style="width:550px"> \
                    <tr> \
                        <td colspan=2> Attracting cycle in the complex $space_name: </td> \
                    </tr> \
                    <tr>\
                        <td style="width:20%;"> Period </td>\
                        <td style="width:80%;"> $(attractor.period) </td>\
                    </tr>\
                    <tr>\
                        <td style="width:20%;"> Multiplier </td>\
                        <td style="width:80%;"> $(attractor.multiplier) </td>\
                    </tr>\
                    <tr>\
                        <td style="width:20%;"> Degree </td>\
                        <td style="width:80%;"> $(attractor.power) </td>\
                    </tr>\
                    <tr>
                        <td style="width:20%;"> Color </td>\
                        <td style="width:80%; color:#$color;"> ██ $(attractor.color) </td>\
                    </tr>\
                    <tr><td style="width:20%;"> Cycle  </td> <td style="width:80%;"> """
        )
    end
    show(io, attractor)
    show(io, m, html"</td> </tr> </table>")
end

function Base.convert(::Type{Attractor{ComplexF64}}, attractor::Attractor{Point})
    cycle = convert(Vector{ComplexF64}, attractor.cycle)
    return Attractor(
        cycle,
        attractor.period,
        attractor.multiplier,
        attractor.power,
        attractor.color,
        attractor.palette,
        attractor.continuous_coloring,
    )
end

const MaybeAttractor = Union{Nothing,Attractor{ComplexF64},Attractor{Point}}

const empty_attractor = nothing

distance(cycle1::Vector{T}, cycle2::Vector{T}) where {T<:PointLike} =
    minimum([distance(z, w) for z in cycle1, w in cycle2])

function is_nearby(z::T, w::T, ε::Real) where {T<:PointLike}
    d = distance(z, w)
    return d < ε ? (true, d, 1) : (false, ε, 0)
end

function is_nearby(z::T, cycle::Vector{T}, ε::Real) where {T<:PointLike}
    for (i, w) in enumerate(cycle)
        d = distance(z, w)
        if d < ε
            return true, d, i
        end
    end
    return false, ε, 0
end

function convergence_time(
    f::Function,
    z::T,
    c::T,
    attractors::Vector{Attractor{T}},
    ε::Float64,
    max_iterations::Int,
) where {T<:PointLike}
    for iteration in 0:max_iterations
        for attractor in attractors
            near, d, shift = is_nearby(z, attractor.cycle, ε)
            near && (return iteration, d, attractor, shift)
        end

        z = f(z, c)
    end

    return max_iterations + 1, ε, empty_attractor, 0
end

struct CloseBy
    iterations::Int
    z::PointLike
    w::PointLike
end

const MaybeCloseBy = Union{Nothing,CloseBy}

function period_multiple_apart(f, z, c, ε, max_iterations)
    slow = fast = z

    for n in 1:max_iterations
        slow = f(slow, c)
        fast = f(f(fast, c), c)

        distance(slow, fast) <= ε && return CloseBy(n, slow, fast)
    end

    return nothing
end

function iterate_until_close(f, z, w, c, ε, max_iterations)
    for n in 1:max_iterations
        z = f(z, c)

        distance(z, w) < ε && return CloseBy(n, z, w)
    end

    return nothing
end

function iterate_both_until_close(f, z, w, c, ε, max_iterations)
    for n in 1:max_iterations
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

function multiplier(f, z0, c, _, ε, max_iterations)
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

function attractor_index(z, julia, d_system, options)
    coloring_data = julia.coloring_data
    attractors = coloring_data.attractors

    f = d_system.map
    N = options.max_iterations
    ε = options.convergence_radius
    c = julia.parameter

    attractor_index = 0
    iterations = 0

    z = convert(attractor_type(coloring_data), z)

    while attractor_index == 0 && iterations <= N
        for (i, attractor) in enumerate(attractors)
            near, _, _ = is_nearby(z, attractor.cycle, ε)
            near && (attractor_index = i)
        end

        z = f(z, c)
        iterations += 1
    end

    return attractor_index
end

# --------------------------------------------------------------------------------------- #
# Finding Attractors
# --------------------------------------------------------------------------------------- #

function coeffs(polynomial, z)
    polynomial = expand(polynomial)
    d = Symbolics.degree(polynomial)

    coefficients::Vector{ComplexF64} = [Symbolics.coeff(polynomial, z^d) for d in 1:d]
    prepend!(coefficients, Symbolics.substitute(polynomial, Dict(z => 0)))
    return convert.(ComplexF64, coefficients)
end

function critical_points(func, parameter)
    # We use Nemo.jl only to extract numerators and denominators of rational functions.
    # Using Symbolics.jl seemed to require checking some edge cases.

    # Create the Complex Rational Funtions Field
    CC = Nemo.ComplexField()
    T, z = Nemo.polynomial_ring(CC, "z")
    F = Nemo.fraction_field(T)

    # Transform function into a symbolic expression
    f = func(F(z), CC(parameter))

    # Get coefficients of numerator and denominator
    num_coeffs = convert.(ComplexF64, collect(Nemo.coefficients(Nemo.numerator(f))))
    den_coeffs = convert.(ComplexF64, collect(Nemo.coefficients(Nemo.denominator(f))))

    num = Polynomial(num_coeffs)
    den = Polynomial(den_coeffs)
    d = max(Polynomials.degree(num), Polynomials.degree(den))

    derivative_num = Polynomials.derivative(num) * den - num * Polynomials.derivative(den)
    points = convert.(ComplexF64, roots(derivative_num))

    while length(points) < 2 * d - 2
        pushfirst!(points, ∞)
    end

    return points
end

function attracting_cycle(f, z::T, c, ε, max_iterations) where {T<:PointLike}
    points = period_multiple_apart(f, z, c, ε, max_iterations)
    isnothing(points) && return T[]

    orbiter = points.z
    reference = points.w

    orbit = [orbiter]
    for _ in 1:max_iterations
        orbiter = f(orbiter, c)
        distance(orbiter, reference) <= ε && return orbit
        push!(orbit, orbiter)
    end

    return orbit
end

function get_attractor(f::Function, z::Number; projective::Bool = false, ε::Real = 1e-4)
    hasmethod(f, ComplexF64) || throw("If it is a family of functions, input a parameter.")

    h = extend_family((z, c) -> f(z))

    return get_attractor(h, 0.0im, z; projective = projective, ε = ε)
end

function get_attractor(
    f::Function,
    c::Number,
    z::Number;
    projective::Bool = false,
    ε::Real = 1e-4,
)
    # If it is not a family ignores the parameter
    hasmethod(f, Tuple{ComplexF64,ComplexF64}) ||
        return get_attractors(f, projective = projective, ε = ε)

    c = convert(ComplexF64, c)
    z = convert(ComplexF64, z)

    # Extend the function, in case it was not extended previously
    g = extend_family(f)

    # Iterate critical points to find all attracting cycles
    crit_pts = critical_points(g, c)
    unique_crit_pts = [(z, count(w -> isapprox(z, w), crit_pts)) for z in unique(crit_pts)]

    pt = convert(Point, z)
    orbit = attracting_cycle(g, pt, c, ε / 2, 1000)
    power = 1

    for (crit_z, multiplicity) in unique_crit_pts
        pt = convert(Point, crit_z)
        power *= is_nearby(pt, orbit, ε / 2)[1] ? (multiplicity + 1) : 1
    end

    plane_orbit = convert(Vector{ComplexF64}, orbit)
    multiplier = abs(g(plane_orbit, c))

    projective || (orbit = plane_orbit)

    return Attractor(orbit, multiplier, power, 1, 1)
end

"""
    get_attractors(f::Function; projective::Bool = false, ε::Real = 1e-4)
    get_attractors(f::Function, c::Number; projective::Bool = false, ε::Real = 1e-4)
    get_attractors(viewer::Viewer)

Find the attracting cycles of a complex map or family. If the input is a `Viewer`, it will \
use the map/family and parameter currently shown in the `viewer`.
"""
function get_attractors(f::Function; projective::Bool = false, ε::Real = 1e-4)
    hasmethod(f, ComplexF64) || throw("If it is a family of functions, input a parameter.")

    h = extend_family((z, c) -> f(z))

    return get_attractors(h, 0.0im; projective = projective, ε = ε)
end

function get_attractors(f::Function, c::Number; projective::Bool = false, ε::Real = 1e-4)
    # If it is not a family ignores the parameter
    hasmethod(f, Tuple{ComplexF64,ComplexF64}) ||
        return get_attractors(f, projective = projective, ε = ε)

    c = convert(ComplexF64, c)

    # Extend the function, in case it was not extended previously
    g = extend_family(f)

    # Iterate critical points to find all attracting cycles
    crit_pts = critical_points(g, c)
    unique_crit_pts = [(z, count(w -> isapprox(z, w), crit_pts)) for z in unique(crit_pts)]

    orbits = Vector{Point}[]
    for (z, multiplicity) in unique_crit_pts
        pt = convert(Point, z)
        orbit = attracting_cycle(g, pt, c, ε / 2, 1000)

        if length(orbit) > 0
            is_repeat = any([distance(orbit, o) < ε for o in orbits])
            !is_repeat && push!(orbits, orbit)
        end
    end

    powers = ones(Int, length(orbits))
    for (i, orbit) in enumerate(orbits)
        for (z, multiplicity) in unique_crit_pts
            pt = convert(Point, z)
            powers[i] *= is_nearby(pt, orbit, ε / 2)[1] ? (multiplicity + 1) : 1
        end
    end

    plane_orbits = convert.(Vector{ComplexF64}, orbits)
    multipliers = [abs(g(o, c)) for o in plane_orbits]

    projective || (orbits = plane_orbits)

    n = length(orbits)
    a = [
        Attractor(o, m, p, i, n) for
        (i, (o, m, p)) in enumerate(zip(orbits, multipliers, powers))
    ]

    return a
end

function rays(func, parameter, show_rays, options)
    if show_rays == false
        return []
    end
    @variables z, c
    f = func(z, c) |> Symbolics.value
    parameter = convert(ComplexF64, parameter)
    coefficients = coeffs(substitute(f, Dict(c => parameter)), z)

    if show_rays == :auto
        periods = [options.period]
        return Rays.auto_rays(coefficients, periods, options.pullbacks)
    elseif show_rays == :all
        return Rays.all_periodic(coefficients, options)
    else
        return Rays.rays(coefficients, show_rays)
    end
end

# --------------------------------------------------------------------------------------- #
# Coloring Algorithms
# --------------------------------------------------------------------------------------- #

const DEFAULT_COLOR_LEVEL = 0.5
const DEFAULT_COLOR = twilight_RGB[255]

to_color(::Nothing) = DEFAULT_COLOR

function to_color(approach::OrbitData)
    depth = mod(approach.preperiod / approach.period / 64.0, 1.0)
    return twilight_RGB[round(Int, 509 * depth + 1)]
end

to_color(::Integer, ::Number, ::Nothing, ::Integer, ::Val{T}) where {T} = DEFAULT_COLOR

function to_color(
    iterations::Integer,
    d::Number,
    attractor::Attractor,
    shift::Integer,
    ::Val{:depth},
)
    d = max(d, 0)
    depth = iterations / attractor.period

    if attractor.continuous_coloring
        if d > 0
            if attractor.multiplier < 1
                depth -= log(attractor.multiplier, d)
            end
            if attractor.power > 1
                depth -= log(attractor.power, -log(d))
            end
        end

        return attractor.palette[round(Int, 509 * mod(depth / 64.0, 1.0) + 1)]
    end

    return attractor.palette[127*mod(round(Int, depth), 4)+95]
end

function to_color(
    iterations::Integer,
    d::Number,
    attractor::Attractor,
    shift::Integer,
    ::Val{:mod},
)
    d = max(d, 0)

    spacing = 255 / attractor.period
    depth = mod(iterations - shift, attractor.period) * spacing + 1

    return attractor.palette[floor(Int, depth)]
end

function escape_time(f, z, c, a, ε, N)
    times = [convergence_time(f, zi, c, a, ε, N) for zi in z]
    ultimate_index = findmin(map(t -> t[1], times))[2]
    ultimate_escape = times[ultimate_index]
    # return max_iterations + 1, ε, empty_attractor, 0
    return to_color(ultimate_escape..., Val(:depth))
end
mod_period(f, z, c, a, ε, N) = to_color(convergence_time(f, z, c, a, ε, N)..., Val(:mod))
convergence_color(f, z, c, a, ε, N) = to_color(multiplier(f, z, c, a, ε, N))

struct ColoringData{T}
    method::Function
    attractors::Vector{Attractor{T}}
    update_attractors::Bool
end

attractor_type(::ColoringData{T}) where {T} = T
attractor_type(is_projective::Bool) = is_projective ? Point : ComplexF64

const base_hue_chroma = [
    #RGB 5794D0 87 148 208 blue
    (36.822489073365745, 266.89189631383687),
    #RGB CA736C 202 115 108 red
    (38.35791098655309, 29.49442824523494),
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
        return :twilight, twilight_RGB
    elseif n_colors < 9
        chroma, hue = base_hue_chroma[color_index]
    else
        chroma, hue = 40, mod(360 * (color_index - 1) / n + 30, 360)
    end

    gradient = [LCHab(20 * cospi(t) + 50, chroma, hue) for t in range(0.0, 2.0, 510)]
    return LCHab(50, chroma, hue), convert.(RGBA{Float64}, gradient)
end

function generate_palette!(attractor::Attractor{T}, color) where {T<:PointLike}
    color = convert.(LCHab, color)
    chroma, hue = color.c, color.h

    palette = [LCHab(20 * cospi(t) + 50, chroma, hue) for t in range(0.0, 2.0, 510)]
    palette = convert.(RGBA{Float64}, palette)

    attractor.color = color
    attractor.palette = palette

    return attractor
end

# --------------------------------------------------------------------------------------- #
# Views of Mandelbrot and Julia Sets
# --------------------------------------------------------------------------------------- #

abstract type View end

mutable struct ColoringScheme
    color::Union{Symbol,LCHab}
    palette::Vector{RGBA{Float64}}
    continuous_coloring::Bool
end

mutable struct Options
    convergence_radius::Float64
    max_iterations::Int
    orbit_length::Int
    critical_length::Int
    compact_view::Bool
    is_family::Bool
    projective_metrics::Tuple{Bool,Bool}
    coloring_methods::Tuple{Symbol,Symbol}
    coloring_schemes::Vector{ColoringScheme}
    pullbacks::Int
    period::Int
    drag_setting::Symbol
end

mutable struct MandelView <: View
    center::ComplexF64
    diameter::Float64
    pixels::Int

    init_center::ComplexF64
    init_diameter::Float64

    colors::Observable{Matrix{RGBA{Float64}}}
    points::Observable{Vector{ComplexF64}}
    marks::Vector{Observable{Vector{ComplexF64}}}
    rays::Vector{Observable{Vector{ComplexF64}}}
    line_refs::Vector{Any}
    mark_line_refs::Vector{Any}
    refresh_rays::Function
    refresh_marks::Function

    coloring_data::ColoringData
    show_rays::Any

    function MandelView(; center, diameter, pixels, coloring_data)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        colors = zeros(RGBA{Float64}, pixels, pixels)
        points = ComplexF64[center]
        marks = Observable{Vector{ComplexF64}}[]
        rays = Observable{Vector{ComplexF64}}[]

        return new(
            center,
            diameter,
            pixels,
            center,
            diameter,
            colors,
            points,
            marks,
            rays,
            Vector{Any}[],
            Vector{Any}[],
            () -> nothing,
            () -> nothing,
            coloring_data,
            0,
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

    colors::Observable{Matrix{RGBA{Float64}}}
    points::Observable{Vector{ComplexF64}}
    marks::Vector{Observable{Vector{ComplexF64}}}
    rays::Vector{Observable{Vector{ComplexF64}}}
    line_refs::Vector{Any}
    mark_line_refs::Vector{Any}
    refresh_rays::Function
    refresh_marks::Function
    coloring_data::ColoringData
    show_rays::Any

    function JuliaView(; center, diameter, parameter, pixels, coloring_data, show_rays)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        colors = zeros(RGBA{Float64}, pixels, pixels)
        points = ComplexF64[center]
        marks = []
        rays = Observable{Vector{ComplexF64}}[]

        return new(
            center,
            diameter,
            parameter,
            pixels,
            center,
            diameter,
            colors,
            points,
            marks,
            rays,
            Vector{Any}[],
            Vector{Any}[],
            () -> nothing,
            () -> nothing,
            coloring_data,
            show_rays,
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
    @inbounds for i in 1:pxs
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

    @inbounds for j in 1:view.pixels
        futures[j] = Threads.@spawn mandel_slice!(
            view.colors[],
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
    @inbounds for i in 1:pxs
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

    @inbounds for j in 1:view.pixels
        futures[j] = Threads.@spawn julia_slice!(
            view.colors[],
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
    notify(view.colors)
    notify(view.points)
    for marks in view.marks
        notify(marks)
    end

    for ray in view.rays
        notify(ray)
    end

    return view
end

function set_marks!(d_system::DynamicalSystem, julia::JuliaView, options::Options)
    c = julia.parameter
    cpts = [z for z in d_system.critical_point(c)]
    for (marks, crpt) in zip(julia.marks, cpts)
        marks[] = orbit(d_system.map, crpt, julia.parameter, options.critical_length - 1)
    end
end

function pick_parameter!(
    julia::JuliaView,
    mandel::MandelView,
    d_system::DynamicalSystem,
    options::Options,
    point,
)
    julia.parameter = point
    set_marks!(d_system, julia, options)
    julia.refresh_marks()

    julia.marks[] = orbit(
        d_system.map,
        d_system.critical_point(julia.parameter),
        julia.parameter,
        options.critical_length - 1,
    )
    julia.rays = []
    julia.refresh_rays()

    mandel.points[] = [julia.parameter]

    coloring_data = julia.coloring_data
    if coloring_data.update_attractors
        T = attractor_type(coloring_data)
        attractor_list = get_attractors(
            d_system.map,
            julia.parameter,
            projective = (T == Point),
            ε = options.convergence_radius / 1000,
        )
        julia.coloring_data = ColoringData{T}(
            coloring_data.method,
            attractor_list,
            coloring_data.update_attractors,
        )
        sync_schemes!(options, julia.coloring_data.attractors)
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
    else
        left_axis.width = nothing
        left_axis.height = nothing
        left_axis.valign = :center
        left_axis.halign = :center

        figure[1, 1][1, 1] = left_axis
        figure[1, 1][1, 2] = right_axis
    end

    for axis in [left_axis, right_axis]
        hidedecorations!(axis)
        deregister_interaction!(axis, :rectanglezoom)
        deregister_interaction!(axis, :limitreset)
        deregister_interaction!(axis, :dragpan)
    end

    return Frame(left_axis, mandel, Dict()), Frame(right_axis, julia, Dict())
end

function delete_plots!(frame::Frame)
    empty!(frame.axis)
    empty!(frame.events)

    view = frame.view[]

    # Clear old listeners (point_vectors, mark_vectors)
    empty!(view.points.listeners)
    # empty!(view.marks.listeners) 
    # TODO: marks and rays
end

function create_plot!(frame::Frame)
    delete_plots!(frame)
    view = frame.view[]

    heatmap!(
        frame.axis,
        view.colors,
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

    view.mark_line_refs = []
    function draw_crit_orbits()
        for r in view.mark_line_refs
            delete!(frame.axis, r)
        end
        for marks in view.marks
            mark_vectors = lift(marks) do zs
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
        end
    end
    view.refresh_marks = draw_crit_orbits
    draw_crit_orbits()

    view.line_refs = []
    function rays_callback()
        for r in view.line_refs
            delete!(frame.axis, r)
        end
        view.line_refs = []
        for ray in view.rays
            ray_vectors = lift(ray) do zs
                xs, ys = to_pixel_space(view, zs)
                return Point2f.(xs, ys)
            end

            push!(
                view.line_refs,
                lines!(
                    frame.axis,
                    ray_vectors,
                    color = (:yellow, 0.5),
                    inspectable = false,
                ),
            )
        end
    end

    rays_callback()
    view.refresh_rays = rays_callback

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

    heatmap!(ax, view.colors[], colormap = :twilight, colorrange = (0.0, 1.0))

    Makie.save(filename, fig)
end

function add_frame_events!(
    figure::Figure,
    frame::Frame,
    topframe::Frame,
    d_system::DynamicalSystem,
    options::Options,
    julia::JuliaView,
)
    axis = frame.axis
    scene = axis.scene

    # Mouse Events
    drag_mode = :notdragging
    dragstart = Point2f(0.0)
    dragend = Point2f(0.0)

    to_world_at_start = z -> to_world(scene, z)

    is_zooming = false
    zooming = Timer(_ -> zoom!(frame, d_system, options), 0.1)

    is_topframe = frame == topframe
    z_level = is_topframe ? 10 : 0

    function set_red_point_useing_mouse_position()
        mp = mouseposition_px(scene)
        view = frame.view[]
        point = to_complex_plane(view, to_world_at_start(mp))
        if view isa MandelView
            pick_parameter!(julia, view, d_system, options, point)
        elseif view isa JuliaView
            pick_orbit!(view, d_system, options, point)
        end
    end

    on(events(scene).mousebutton) do event
        is_zooming && return Consume(false)

        mp = mouseposition_px(scene)
        view = frame.view[]

        if event.button == Mouse.right
            if event.action == Mouse.press &&
               is_mouseinside(axis) &&
               (is_topframe || !is_mouseinside(topframe.axis))
                drag_mode = :rightclick
                to_world_at_start = z -> to_world(scene, z)
                dragstart = to_world_at_start(mp)

            elseif event.action == Mouse.release && (drag_mode == :rightclick)
                dragend = to_world_at_start(mp)
                view.center +=
                    to_complex_plane(view, dragstart) - to_complex_plane(view, dragend)
                translate!(scene, 0, 0, z_level)
                update_view!(view, d_system, options)
                reset_limits!(axis)
                drag_mode = :notdragging
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

                elseif ispressed(scene, Keyboard.left_shift | Keyboard.right_shift) &&
                       view == julia &&
                       options.coloring_methods[2] != :preperiod
                    z = to_complex_plane(julia, mouseposition(scene))
                    i = attractor_index(z, julia, d_system, options)

                    i == 0 && (return Consume(true))

                    change_color!(figure, julia, i, d_system, options)
                    return Consume(true)
                elseif options.drag_setting == :both ||
                       options.drag_setting == :dynamic_only && view isa JuliaView
                    set_red_point_useing_mouse_position()
                    drag_mode = :leftclick
                else
                    set_red_point_useing_mouse_position()
                end
            elseif event.action == Mouse.release && drag_mode == :leftclick
                set_red_point_useing_mouse_position()
                drag_mode = :notdragging
            end
        end

        return Consume(false)
    end

    on(events(scene).mouseposition) do event
        if (drag_mode == :rightclick)
            # TODO: figure out differences between mouseposition mouseposition_px because we use both in a way that I believe causes issues.
            mp = mouseposition(scene)
            translate!(scene, mp - dragstart..., z_level)
        elseif (drag_mode == :leftclick)
            set_red_point_useing_mouse_position()
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

        if ispressed(scene, (Keyboard.left_control | Keyboard.right_control) & Keyboard.p)
            if is_mouseinside(axis) &&
               !is_zooming &&
               (is_topframe || !is_mouseinside(topframe.axis))
                println(view.points.val[1])
                return Consume(true)
            end
        end

        if ispressed(scene, Keyboard.c) &&
           view == julia &&
           options.coloring_methods[2] != :preperiod
            if is_mouseinside(axis) &&
               !is_zooming &&
               (is_topframe || !is_mouseinside(topframe.axis))
                z = to_complex_plane(julia, mouseposition(scene))
                i = attractor_index(z, julia, d_system, options)

                i == 0 && (return Consume(true))

                julia.coloring_data.attractors[i].continuous_coloring =
                    !julia.coloring_data.attractors[i].continuous_coloring

                options.coloring_schemes[i].continuous_coloring =
                    julia.coloring_data.attractors[i].continuous_coloring

                update_view!(julia, d_system, options)
                return Consume(true)
            end
        end
    end
end

function add_buttons!(
    figure,
    left_frame,
    right_frame,
    mandel,
    julia,
    d_system,
    options,
    show_rays,
)
    layout = GridLayout(figure[2, 1], tellwidth = false)
    button_shift = options.is_family ? 2 : 0

    labels = Dict(
        :max_iter => Label(layout[1, button_shift+1], "Maximum\nIterations:"),
        :orbit_len => Label(layout[1, button_shift+3], "Orbit\nLength:"),
        :critical_length => Label(layout[1, button_shift+5], "Critical\nOrbit Length:"),
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

    if show_rays == :auto || show_rays == :all
        labels[:pullbacks] = Label(layout[1, button_shift+9], "Pullbacks:")
        inputs[:pullbacks] = Textbox(
            layout[1, button_shift+10],
            width = 60,
            placeholder = string(options.pullbacks),
            validator = Int,
        )
        labels[:period] = Label(layout[1, button_shift+11], "Max\n Period:")
        inputs[:period] = Textbox(
            layout[1, button_shift+12],
            width = 60,
            placeholder = string(options.period),
            validator = Int,
        )
        inputs[:rays] = Button(layout[1, button_shift+13], label = "R", halign = :left)

        on(inputs[:pullbacks].stored_string) do s
            options.pullbacks = parse(Int, s)
        end
        on(inputs[:period].stored_string) do s
            options.period = parse(Int, s)
        end

        on(inputs[:rays].clicks, priority = 200) do event
            inputs[:rays].label = "X"

            new_rays = rays(d_system.map, julia.parameter, julia.show_rays, options)
            julia.rays = [Observable(ray) for ray in new_rays]
            julia.refresh_rays()

            inputs[:rays].label = "R"
        end
    end

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

        set_marks!(d_system, julia, options)
    end

    on(inputs[:convergence_radius].stored_string) do s
        options.convergence_radius = parse(Float64, s)

        coloring_data = julia.coloring_data
        if coloring_data.update_attractors
            T = attractor_type(coloring_data)
            attractor_list = get_attractors(
                d_system.map,
                julia.parameter,
                projective = (T == Point),
                ε = options.convergence_radius / 1000,
            )
            julia.coloring_data = ColoringData{T}(
                coloring_data.method,
                attractor_list,
                coloring_data.update_attractors,
            )
            sync_schemes!(options, julia.coloring_data.attractors)
        end

        update_view!(julia, d_system, options)
        options.is_family && update_view!(mandel, d_system, options)
    end

    if options.is_family
        inputs[:switch_layout] = Button(layout[1, 1], label = "↰", halign = :left)
        inputs[:switch_positions] = Button(layout[1, 2], label = "↔", halign = :left)

        on(inputs[:switch_layout].clicks, priority = 200) do event
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

        on(inputs[:switch_positions].clicks, priority = 200) do event
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
    options = [:escape_time, :convergence_time, :mod_period, :preperiod]
    mandel in options || throw("Invalid Mandelbrot set coloring method")
    julia in options || throw("Invalid Julia set coloring method")

    # Mandelbrot can't use attractors so change to default
    mandel == :convergence_time && (mandel = :preperiod)
    mandel == :mod_period && (mandel = :preperiod)
    return mandel, julia
end

function get_coloring_data(map, c, coloring_method, projective_metric)
    if coloring_method == :escape_time
        method = escape_time
        a = get_attractor(map, c + 100, ∞, projective = projective_metric)
        attractors = [a]
        update_attractors = false
    elseif coloring_method == :convergence_time
        method = escape_time
        attractors = get_attractors(map, c, projective = projective_metric)
        update_attractors = true
    elseif coloring_method == :mod_period
        method = mod_period
        attractors = get_attractors(map, c, projective = projective_metric)
        update_attractors = true
    elseif coloring_method == :preperiod
        method = convergence_color
        attractors = Attractor{attractor_type(projective_metric)}[]
        update_attractors = false
    else
        throw("Invalid `coloring_method`")
    end

    return ColoringData(method, attractors, update_attractors)
end

function get_scheme(a::Attractor{T}) where {T}
    return ColoringScheme(a.color, a.palette, a.continuous_coloring)
end

function store_schemes!(options::Options, attractors::Vector{Attractor{T}}) where {T}
    options.coloring_schemes = [get_scheme(a) for a in attractors]
    return options
end

function sync_schemes!(options::Options, attractors::Vector{Attractor{T}}) where {T}
    schemes = options.coloring_schemes

    for (s, a) in zip(schemes, attractors)
        a.color = s.color
        a.palette = s.palette
        a.continuous_coloring = s.continuous_coloring
    end

    start = length(schemes) + 1

    for attractor in attractors[start:end]
        push!(schemes, get_scheme(attractor))
    end

    return options
end

"""
    Viewer(f; <keyword arguments>)

Create a `Viewer` that plots the Mandelbrot and the Julia sets associated with the \
function `f`. `f` can have one or two inputs, the second being the parameter. By default, \
the `Viewer` opens as a separate window.

# Examples

```julia-repl
Viewer((z, c) -> z^2 + c, mandel_center = -0.5)
```

```julia-repl
f(z, λ) = z^2 + λ / z^2
crit(λ) = λ^(1 / 4)
Viewer(f; crit = crit, mandel_diameter = 1.0)
```

# Shortcuts

  - `Right Click Drag`: Pans view.
  - `Mouse Scroll`: Zooms in or out.
  - `Left Click` (parameter space): Chooses parameter.
  - `Left Click` (dynamical space): Chooses orbit initial point.
  - `C` (dynamical space): Change between discrete and continuous coloring.
  - `Shift + Left Click`: Change color of the basin of attraction.
  - `Ctrl + Left Click`: Resets view to initial `center` and `diameter`.
  - `Ctrl + S`: Saves view that is under the mouse pointer (files in `./imgs`).
  - `Ctrl + P`: Prints the value of the red point in the view under the pointer.

# Arguments

  - `crit = 0.0im`: Function that gives a critical point for each parameter. Used to plot \
    the Mandelbrot set. If it is a constant function, you can just input the constant \
    directly.
  - `c = 0.0im`: Initial parameter used to plot the Julia set.
  - `mandel_center = 0.0im`: Initial center of the Mandelbrot plot.
  - `mandel_diameter = 4.0`: Initial diameter of the Mandelbrot plot.
  - `julia_center = 0.0im`: Initial center of the Julia plot.
  - `julia_diameter = 4.0`: Initial diameter of the Julia plot.
  - `grid_width = 800`: Width (and height) of the grid of complex numbers used to plot \
    sets.
  - `compact_view = false`: If 'true' one of the plots is show as an inset plot, if \
    `false` they are shown side-by-side.
  - `show_rays = false`: Rays can only be computed for polynomials. Only the dynamic \
    rays can be computed as yet. If `false`, no  rays are shown. If a vector of \
    Rational64 is given, then the orbits of those rays are displayed. If `:all` \
    is given then a button will be added to compute all the rays up to a period and \
    pullbacks. If `:auto` is given, then those rays are then filtered by whether they \
    land at a cut point, and a lamination is printed.
  - `left_click_drag = :dynamic_only`: By default, in the dynamic plain, the red \
    point will be continuously updated to the mouse postion if the left button \
    remains pressed. This may be undesiorable for performace reasons. Set to `:neither` \
    to disable. Alternately, it may be tolerable to enable this in both plains with \
    `:both`.

# Coloring Method Options

For the options below, if you want to set different values for the Mandelbrot and Julia \
set views, set the option to be a tuple with the respective values.

  - `coloring_method = :escape_time`: Chooses the coloring method for both plots. \
    The options are `:escape_time`, `:convergence_time`, `:mod_period`, `:preperiod`. \
    More details below.
  - `projective_metric = false`: Determines which metric will be used to determine \
    distances, the complex plane metric, or the metric on the projective line. The \
    distance between ∞ and a finite point in the plane metric is the inverse of its absolute value.

# Coloring Methods

  - `:escape_time` (default): colors according how fast the point approaches ∞;
  - `:convergence_time`: computes all attracting cycles in advance, and then computes how \
    fast each point converges to one of those attractors (uses different color gradients \
    for each attractor). This option is only available for the Julia set, and it will \
    default to `:preperiod` in the Mandelbrot set case.
  - `:mod_period`: similar to `:convergence_time` but instead of using the `preperiod` \
    (i.e. how fast it converges to attractor) as the "color depth" it uses \
    `mod(preperiod, period)`.
  - `:preperiod`: finds the attracting cycle for each point separately using \
    Floyd's cycle-finding algorithm and uses `preperiod / period` as "color depth";
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
        coloring_method = :escape_time,
        projective_metric = false,
        show_rays = false,
        left_click_drag = :dynamic_only,
    )
        # Put options in standard form
        projective_metrics = make_tuple(projective_metric)
        coloring_methods = fix_criteria(make_tuple(coloring_method))

        # Check if it is a family of maps
        is_family = hasmethod(f, Tuple{ComplexF64,ComplexF64})

        # Create Viewer Data
        d_system = DynamicalSystem(f, crit)
        options = Options(
            1e-3,
            200,
            1,
            1,
            compact_view,
            is_family,
            projective_metrics,
            coloring_methods,
            ColoringScheme[],
            0,
            1,
            left_click_drag,
        )
        figure = Figure(size = (800, 850))

        mandel_coloring =
            get_coloring_data(d_system.map, c, coloring_methods[1], projective_metrics[1])
        julia_coloring =
            get_coloring_data(d_system.map, c, coloring_methods[2], projective_metrics[2])

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
            show_rays = show_rays,
        )

        test_crits = d_system.critical_point(c)
        test_crits = test_crits isa Number ? (test_crits,) : test_crits
        julia.marks = [Observable(ComplexF64[]) for _ in test_crits]

        store_schemes!(options, julia_coloring.attractors)

        set_marks!(d_system, julia, options)

        julia.rays = []
        pick_orbit!(julia, d_system, options, julia.points[][begin])

        left_frame, right_frame = create_frames!(figure, options, mandel, julia)

        if is_family
            add_frame_events!(figure, right_frame, left_frame, d_system, options, julia)
            add_frame_events!(figure, left_frame, left_frame, d_system, options, julia)

            create_plot!(left_frame)
            update_view!(mandel, d_system, options)
            mandel.points[] = [julia.parameter]
        else
            add_frame_events!(figure, right_frame, right_frame, d_system, options, julia)
        end

        create_plot!(right_frame)
        update_view!(julia, d_system, options)

        colgap!(content(figure[1, 1]), 10)
        inputs = add_buttons!(
            figure,
            left_frame,
            right_frame,
            mandel,
            julia,
            d_system,
            options,
            show_rays,
        )

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

function get_attractors(viewer::Viewer)
    attractors =
        viewer.julia.coloring_data.update_attractors ?
        viewer.julia.coloring_data.attractors :
        get_attractors(viewer.d_system.map, viewer.julia.parameter)

    return attractors
end

"""
    get_parameter(viewer::Viewer)

Get the parameter used to plot the Julia set in `viewer`.
"""
get_parameter(viewer::Viewer) = viewer.julia.parameter

function change_color!(julia, i, chroma, hue, d_system, options)
    isempty(julia.coloring_data.attractors) &&
        throw("Julia set plot is not using attractors for coloring")

    color = LCHab{Float64}(50, chroma, hue)
    attractors = julia.coloring_data.attractors
    attractor = generate_palette!(attractors[i], color)

    options.coloring_schemes[i].color = attractor.color
    options.coloring_schemes[i].palette = attractor.palette
    options.coloring_schemes[i].continuous_coloring = attractor.continuous_coloring

    update_view!(julia, d_system, options)
    return nothing
end

function change_color!(figure, julia, i, d_system, options)
    ax = PolarAxis(figure[1, 1], width = Relative(1.01), height = Relative(1.01))

    ax.rzoomlock = true
    ax.thetazoomlock = true
    ax.fixrmin = true
    ax.r_translation_button = false
    ax.theta_translation_button = false

    chromas = 0:100
    hues = 0:360
    cs = [LCHab(50, c, h) for h in hues, c in chromas]

    surface!(
        ax,
        0 .. 2pi,
        0 .. 100,
        zeros(size(cs)) .+ 100,
        color = cs,
        shading = NoShading,
    )
    ax.gridz[] = 200

    attractors = julia.coloring_data.attractors
    colors = [attractor.color for attractor in attractors if !isa(attractor.color, Symbol)]

    if !isempty(colors)
        chromas = [color.c for color in colors]
        hues = deg2rad.([color.h for color in colors])
        scatter!(
            ax,
            hues,
            chromas,
            zeros(size(colors)) .+ 301,
            color = colors,
            markersize = 8,
        )
        scatter!(
            ax,
            hues,
            chromas,
            zeros(size(colors)) .+ 300,
            color = :black,
            markersize = 16,
        )
    end

    tightlimits!(ax)
    hidedecorations!(ax)

    colorpicker = on(events(ax.scene).mousebutton, priority = 300, weak = true) do event
        mp = mouseposition(ax.scene)
        if event.button == Mouse.left && event.action == Mouse.press
            w = complex(mp[1], mp[2])
            r = abs(w)
            θ = rad2deg(angle(w))
            change_color!(julia, i, r, θ, d_system, options)

            delete!(ax)
            off(colorpicker)
        end

        return Consume(true)
    end

    return nothing
end

end
