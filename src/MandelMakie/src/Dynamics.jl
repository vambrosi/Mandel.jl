using StaticArraysCore, LinearAlgebra

const DEFAULT_VALUE = 0.5

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

    for iter in 1:max_iter
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

    for iter in 1:max_iter
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

const Point = MVector{2, ComplexF64}

const Mobius = SMatrix{2, 2, ComplexF64}

function antipodal_hyperbolic(fixed_point, scale_factor)
    pt_normalized = normalize(fixed_point)
    u = pt_normalized[1]
    v = pt_normalized[2]

    return Mobius(
        abs2(u) + scale_factor * abs2(v),
        (1 - scale_factor) * conj(u) * v,
        (1 - scale_factor) * u * conj(v),
        scale_factor * abs2(u) + abs2(v)
    )
end

function distance(pt1::Point, pt2::Point)
	return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

function multiplier(
	f::Function,
	pt::Point,
	c::Point,
	ε::Real,
	max_iter::Integer
)
    tortoise = hare = pt

	iter = 1
	period_multiple = 1
	while iter <= max_iter
		tortoise = f(tortoise, c)
		hare = f(f(hare, c), c)

		if distance(hare, tortoise) <= ε
			period_multiple = iter
			break
		end

		iter += 1
	end

	if iter == max_iter + 1
		return 0.5
	end

	μ = 1
	iter = 1
	period = 1
	while iter <= period_multiple
		tortoise = f(tortoise, c)
		# μ *= df(tortoise, c)

		if distance(hare, tortoise) <= ε
			period = iter
			break
		end

		iter += 1
	end

	if iter == period_multiple + 1
		return 0.5
	end

	tortoise = pt
	iter = 0
	preperiod = 0
	while iter <= max_iter
		if distance(hare, tortoise) <= ε
			preperiod = iter
			break
		end

		tortoise = f(tortoise, c)
		hare = f(hare, c)
		iter += 1
	end

	if iter == max_iter + 1
		return 0.5
	end

	return mod((preperiod / period) / 64.0, 1.0)
end