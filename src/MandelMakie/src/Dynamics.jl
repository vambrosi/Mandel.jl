using StaticArraysCore, LinearAlgebra

const DEFAULT_VALUE = 0.5

# --------------------------------------------------------------------------------------- #
# Complex Plane Definitions
# --------------------------------------------------------------------------------------- #

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

"""
    Point(z::ComplexF64, w::ComplexF64)

A `Point` in the complex line represented by a pair of complex numbers `z` and `w`.
"""
const Point = MVector{2,ComplexF64}

"""
    Mobius(z11::ComplexF64, z21::ComplexF64, z12::ComplexF64, z22::ComplexF64)

A `Mobius` transformation represented by the 2x2 matrix:
```
z11 z12
z21 z22
```
"""
const Mobius = SMatrix{2,2,ComplexF64}

function to_complex_plane(pt::Point)
    return pt[1] / pt[2]
end

"""
    hyperbolic(fixed_pt1, fixed_pt2, scale_factor)

Mobius transformation with two fixed points which is conjugate to the dilation \
`f(z) = scale_factor * z`. The conjugation maps `(fixed_pt1, fixed_pt2)` to `(0, ∞)`.
"""
function hyperbolic(fixed_pt1, fixed_pt2, scale_factor)
    pt1 = normalize(fixed_pt1)
    pt2 = normalize(fixed_pt2)

    return Mobius(
        pt1[1] * pt2[2] - scale_factor * pt2[1] * pt1[2],
        (1 - scale_factor) * pt1[2] * pt2[2],
        (scale_factor - 1) * pt1[1] * pt2[1],
        scale_factor * pt1[1] * pt2[2] - pt2[1] * pt1[2],
    )
end

"""
    distance(pt1::Point, pt2::Point)

Distance between to points on the complex projective line. This function assumes that \
the representatives have unit norm, i.e. that they are inside the 3-sphere.
"""
function distance(pt1::Point, pt2::Point)
    return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

function orbit(f::Function, pt::Point, c::Point, length::Int)
    points = Vector{Point}(undef, length)
    points[1] = pt

    for i = 2:length
        points[i] = f(points[i-1], c)
    end

    return points
end

# --------------------------------------------------------------------------------------- #
# Coloring Algorithms for the Complex Plane
# --------------------------------------------------------------------------------------- #

function escape_time(
    f::Function,
    df::Function,
    z0::Number,
    c::Number,
    esc_radius::Real,
    max_iter::Integer,
)
    z = z0

    for iter = 1:max_iter
        z = f(z, c)

        if abs(z) > esc_radius
            return mod((iter + 1.0 - log2(log(abs(z)))) / 64.0, 1.0)
        end
    end

    return DEFAULT_VALUE
end

function stop_time(
    f::Function,
    df::Function,
    z0::Number,
    c::Number,
    esc_radius::Real,
    max_iter::Integer,
)
    z = z0
    inv_radius = 1 / 1000 * esc_radius

    for iter = 1:max_iter
        z1 = f(z, c)

        if abs(z1 - z) <= inv_radius
            mult = abs(df(z, c))
            eps = abs(z1 - z)

            return mod((iter + 1.0 - log(eps) / log(mult)) / 64.0, 1.0)
        end

        z = z1
    end

    return DEFAULT_VALUE
end

function escape_preperiod(
    f::Function,
    df::Function,
    z0::Number,
    c::Number,
    esc_radius::Real,
    max_iter::Integer,
)
    inv_radius = 1 / (1000 * esc_radius)

    tortoise = hare = z0

    iter = 1
    period_multiple = 1
    while iter <= max_iter

        if abs2(tortoise) > esc_radius
            return mod((iter + 1.0 - log2(log(abs(tortoise)))) / 64.0, 1.0)
        end

        tortoise = f(tortoise, c)
        hare = f(f(hare, c), c)

        if abs2(hare - tortoise) <= inv_radius
            period_multiple = iter
            break
        end

        iter += 1
    end

    if iter == max_iter + 1
        return DEFAULT_VALUE
    end

    μ = 1
    iter = 1
    period = 1
    while iter <= period_multiple
        tortoise = f(tortoise, c)
        μ *= df(tortoise, c)

        if abs2(hare - tortoise) <= inv_radius
            period = iter
            break
        end

        iter += 1
    end

    if iter == period_multiple + 1
        return DEFAULT_VALUE
    end

    tortoise = z0
    iter = 0
    preperiod = 0
    while iter <= max_iter
        if abs2(hare - tortoise) <= inv_radius
            preperiod = iter
            break
        end

        tortoise = f(tortoise, c)
        hare = f(hare, c)
        iter += 1
    end

    if iter == max_iter + 1
        return DEFAULT_VALUE
    end

    ε = abs(hare - tortoise)

    return mod((preperiod / period + 1 - log(ε) / log(abs(μ))) / 64.0, 1.0)
end

# --------------------------------------------------------------------------------------- #
# Coloring Algorithms for the Complex Projective Line
# --------------------------------------------------------------------------------------- #

struct NearbyPoints
    close::Bool
    iterations::Int
    slow::Point
    fast::Point
end

function find_period_multiple(f::Function, pt::Point, c::Point, ε::Real, max_iter::Integer)
    slow = fast = pt

    for period_multiple = 1:max_iter
        slow = f(slow, c)
        fast = f(f(fast, c), c)

        distance(slow, fast) <= ε && return NearbyPoints(true, period_multiple, slow, fast)
    end

    return NearbyPoints(false, max_iter, slow, fast)
end

function iterate_slow_until_close(
    f::Function,
    slow::Point,
    fast::Point,
    c::Point,
    ε::Real,
    max_iter::Integer,
)
    for period = 1:max_iter
        slow = f(slow, c)

        distance(slow, fast) <= ε && return NearbyPoints(true, period, slow, fast)
    end

    return NearbyPoints(false, max_iter, slow, fast)
end

function orbit_until_close(
    f::Function,
    orbiter::Point,
    reference::Point,
    c::Point,
    ε::Real,
    max_iter::Integer,
)
    orbit = [orbiter]
    for _ = 1:max_iter
        orbiter = f(orbiter, c)
        distance(orbiter, reference) <= ε && return orbit
        push!(orbit, orbiter)
    end

    return orbit
end

function iterate_both_until_close(
    f::Function,
    pt1::Point,
    pt2::Point,
    c::Point,
    ε::Real,
    max_iter::Integer,
)
    for iteration = 0:max_iter
        distance(pt1, pt2) <= ε && return NearbyPoints(true, iteration, pt1, pt2)

        pt1 = f(pt1, c)
        pt2 = f(pt2, c)
    end

    return NearbyPoints(false, max_iter, pt1, pt2)
end

function multiplier(f::Function, pt::Point, c::Point, ε::Real, max_iter::Integer)

    first_approach = find_period_multiple(f, pt, c, ε, max_iter)
    !first_approach.close && return DEFAULT_VALUE

    period_data = iterate_slow_until_close(
        f,
        first_approach.slow,
        first_approach.fast,
        c,
        ε,
        max_iter,
    )
    !period_data.close && return DEFAULT_VALUE

    preperiod_data = iterate_both_until_close(f, pt, period_data.fast, c, ε, max_iter)
    !preperiod_data.close && return DEFAULT_VALUE

    return mod((preperiod_data.iterations / period_data.iterations) / 64.0, 1.0)
end

struct FatouIterationDistance
    preperiod::Int
    period::Int
    component_index::Int
end

function convergence_time(
    f::Function,
    pt::Point,
    c::Point,
    attractors,
    ε::Real,
    max_iter::Integer,
)
    pt = normalize(pt)

    for preperiod = 0:max_iter
        for (component_index, attractor) in enumerate(attractors)
            for limit_pt in attractor
                distance(pt, limit_pt) <= ε && return FatouIterationDistance(
                    preperiod,
                    length(attractor),
                    component_index,
                )
            end
        end

        pt = f(pt, c)
    end

    return FatouIterationDistance(max_iter, 1, 0)
end

function attracting_orbit(f::Function, pt::Point, c::Point, ε::Real, max_iter::Integer)
    first_approach = find_period_multiple(f, pt, c, ε, max_iter)
    !first_approach.close && return Point[]

    orbit = orbit_until_close(f, first_approach.slow, first_approach.fast, c, ε, max_iter)
    return orbit
end
