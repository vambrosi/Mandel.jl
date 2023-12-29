module Projective

mutable struct Point
    u::ComplexF64
    v::ComplexF64
end

"Complex number associated with the point in the projective line."
function to_complex(pt::Point)
    if pt.v != 0
        return pt.u / pt.v
    else
        return Inf + Inf*im
    end
end

function orbit(f::Function, pt::Point, iters::Integer)
    orbit = Array{Point}(undef, iters + 1)

    orbit[1] = pt
    for iter in 2:(iters + 1)
        orbit[iter] = f(orbit[iter - 1])
    end

    return orbit
end

function complex_orbit(f::Function, pt::Point, iters::Integer)
    orbit = orbit(f::Function, pt::Point, iters::Integer)

    return to_complex.(orbit)
end

# Auxiliary functions
norm_sqr(pt::Point) = abs2(pt.u) + abs2(pt.v)

# Distance in the projective plane
dist_sqr(pt1::Point, pt2::Point) = abs2(pt1.u * pt2.v - pt1.v * pt2.u) / (norm_sqr(pt1) * norm_sqr(pt2))
dist(pt1::Point, pt2::Point) = sqrt(dist_sqr(pt1, pt2))

# Computes the distance squared and normalizes points
function d_sqr!(pt1::Point, pt2::Point)
    numerator = abs2(pt1.u * pt2.v - pt1.v * pt2.u)
    norm2_1 = norm_sqr(pt1)
    norm2_2 = norm_sqr(pt2)

    norm1 = sqrt(norm2_1)
    norm2 = sqrt(norm2_2)

    pt1.u /= norm1
    pt1.v /= norm1
    pt2.u /= norm2
    pt2.v /= norm2

    return numerator / (norm2_1 * norm2_2)
end

"Takes two coordinate functions and returns the corresponding map in the projective line."
projectivize(f_u::Function, f_v::Function) = (pt::Point -> Point(f_u(pt.u, pt.v), f_v(pt.u, pt.v)))

"""
Finds the period and preperiod of a point up to a level of precision (given by ε).
Returns max_iter + 1 for period (or preperiod) when it is larger than max_iter.
"""
function dynamics_at(f::Function, pt, max_iter::Integer, ε::Real)
    tortoise = hare = pt
    ε_sqr = ε^2

    # In case condition is not triggered
    period_multiple = max_iter + 1

    # Finds first time z_n and z_2n are close
    for iter in 1:max_iter
        tortoise = f(tortoise)
        hare = f(f(hare))

        if d_sqr!(tortoise, hare) < ε_sqr
            period_multiple = iter
            break
        end
    end

    # Exits if n is bigger than max_iter
    if period_multiple == max_iter + 1
        return period_multiple, period_multiple
    end

    # Finds next time they get close
    period = period_multiple + 1
    for iter in 1:period_multiple
        tortoise = f(tortoise)

        if d_sqr!(tortoise, hare) < ε_sqr
            period = iter
            break
        end
    end

    # Finds preperiod
    tortoise = pt
    preperiod = max_iter + 1
    for iter in 0:max_iter
        if d_sqr!(tortoise, hare) < ε_sqr
            preperiod = iter
            break
        end

        tortoise = f(tortoise)
        hare = f(hare)
    end

    return period, preperiod
end

end