using Polynomials, HypertextLiteral, Crayons, Colors
import Nemo

struct DynamicalSystem
    map::Function
    critical_point::Function

    function DynamicalSystem(f::Function, critical_point::Function)
        result = critical_point(0.0im)
        if !(result isa ComplexF64 || result isa Vector{ComplexF64})
            throw(
                "The critical point function must return a ComplexF64 or a Vector{ComplexF64}.\n" *
                "Got result of type $(typeof(result))",
            )
        end

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
        _ -> [c]
    end

    return DynamicalSystem(f, g)
end

mutable struct Attractor{T}
    cycle::Union{T,Vector{T}}
    period::Int
    multiplier::Float64
    power::Int
    color::Union{Symbol,LCHab}
    palette::Vector{RGB24}
    continuous_coloring::Bool
end

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

function Base.convert(::Type{Attractor{ComplexF64}}, attractor::Attractor{ProjectivePoint})
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

const MaybeAttractor = Union{Nothing,Attractor{ComplexF64},Attractor{ProjectivePoint}}

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
    z::Vector{T},
    c::T,
    attractors::Vector{Attractor{T}},
    ε::Float64,
    max_iterations::Int,
) where {T<:PointLike}
    z = copy(z)

    for iteration in 0:max_iterations
        for attractor in attractors
            for zi in z
                near, d, shift = is_nearby(zi, attractor.cycle, ε)
                near && (return iteration, d, attractor, shift)
            end
        end

        for i in eachindex(z)
            z[i] = f(z[i], c)
        end
    end

    return max_iterations + 1, ε, empty_attractor, 0
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

        for i in eachindex(z)
            z = f(z, c)
        end
    end

    return max_iterations + 1, ε, empty_attractor, 0
end

struct CloseBy
    iterations::Int
    z::PointLike
    w::PointLike
end

const MaybeCloseBy = Union{Nothing,CloseBy}

function period_multiple_apart(f, z::Vector{T}, c, ε, max_iterations) where {T<:PointLike}
    slow = copy(z)
    fast = copy(z)

    for n in 1:max_iterations
        for i in eachindex(slow)
            slow[i] = f(slow[i], c)
            fast[i] = f(f(fast[i], c), c)
        end

        if any(distance(s, f) <= ε for (s, f) in zip(slow, fast))
            close = [distance(s, f) <= ε for (s, f) in zip(slow, fast)]
            index = findall(close)[1]
            return CloseBy(n, slow[index], fast[index]), index
        end
    end

    return nothing
end

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

function multiplier(f, z0::Vector{T}, c, _, ε, max_iterations) where {T<:PointLike}
    result = period_multiple_apart(f, z0, c, ε, max_iterations)
    isnothing(result) && return nothing

    (points, index) = result
    points = iterate_until_close(f, points.z, points.w, c, ε, max_iterations)
    isnothing(points) && return nothing

    period = points.iterations
    points = iterate_until_close(f, z0[index], points.w, c, ε, max_iterations)
    isnothing(points) && return nothing

    preperiod = points.iterations
    return OrbitData(preperiod, period)
end

function multiplier(f, z0::T, c, _, ε, max_iterations) where {T<:PointLike}
    return multiplier(f, [z0], c, nothing, ε, max_iterations)
end

function attractor_index(z, julia, d_system, options)
    coloring_data = julia.coloring_data
    attractors = coloring_data.attractors

    f = d_system.map
    N = options.max_iterations
    ε = options.convergence_radius
    c = get_parameter(julia)

    attractor_index = 0
    iterations = 0

    z = convert(get_attractor_type(coloring_data), z)

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

    pt = convert(ProjectivePoint, z)
    orbit = attracting_cycle(g, pt, c, ε / 2, 1000)
    power = 1

    for (crit_z, multiplicity) in unique_crit_pts
        pt = convert(ProjectivePoint, crit_z)
        power *= is_nearby(pt, orbit, ε / 2)[1] ? (multiplicity + 1) : 1
    end

    plane_orbit = convert(Vector{ComplexF64}, orbit)
    multiplier = abs(g(plane_orbit, c, Val(:diff)))

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

    orbits = Vector{ProjectivePoint}[]
    for (z, multiplicity) in unique_crit_pts
        pt = convert(ProjectivePoint, z)
        orbit = attracting_cycle(g, pt, c, ε / 2, 1000)

        if length(orbit) > 0
            is_repeat = any([distance(orbit, o) < ε for o in orbits])
            !is_repeat && push!(orbits, orbit)
        end
    end

    powers = ones(Int, length(orbits))
    for (i, orbit) in enumerate(orbits)
        for (z, multiplicity) in unique_crit_pts
            pt = convert(ProjectivePoint, z)
            powers[i] *= is_nearby(pt, orbit, ε / 2)[1] ? (multiplicity + 1) : 1
        end
    end

    plane_orbits = convert.(Vector{ComplexF64}, orbits)
    multipliers = [abs(g(o, c, Val(:diff))) for o in plane_orbits]

    projective || (orbits = plane_orbits)

    n = length(orbits)
    a = [
        Attractor(o, m, p, i, n) for
        (i, (o, m, p)) in enumerate(zip(orbits, multipliers, powers))
    ]

    return a
end