# --------------------------------------------------------------------------------------- #
# Parsing Functions
# --------------------------------------------------------------------------------------- #
using Polynomials

struct DynamicalSystem
	f::Function
	df_dz::Function
	df_dc::Function
	crit::Function

	function DynamicalSystem(f, crit)
        @variables z, c
		expr = f(z, c)

		df_dz_expr = expand_derivatives(Differential(z)(expr))
		df_dz = build_function(df_dz_expr, z, c, expression=Val{false})

		df_dc_expr = expand_derivatives(Differential(c)(expr))
		df_dc = build_function(df_dc_expr, z, c, expression=Val{false})

        if crit isa Function
            crit_function = crit
        else
            critical_point = convert(ComplexF64, crit)
            crit_function = (_) -> critical_point
        end

		return new(f, df_dz, df_dc, crit_function)
	end
end

function derivative(f)
	@variables z, c
	df = expand_derivatives(Differential(z)(f(z, c)))
	return build_function(df, z, c, expression=Val{false})
end

function to_point_family(f)
	@variables u, v, c, d
	frac = simplify(f(u / v, c / d))
	value = Symbolics.value(frac)

	num, den = Symbolics.arguments(value)
	fu = build_function(num, u, v, c, d, expression=Val{false})
	fv = build_function(den, u, v, c, d, expression=Val{false})

	return (pt_z, pt_c) -> normalize!(Point(fu(pt_z..., pt_c...), fv(pt_z..., pt_c...)))
end

function to_point_map(f)
	@variables c, d
	frac = simplify(f(c / d))
	value = Symbolics.value(frac)

	num, den = Symbolics.arguments(value)
	fu = build_function(num, c, d, expression=Val{false})
	fv = build_function(den, c, d, expression=Val{false})

	return (pt) -> normalize!(Point(fu(pt...), fv(pt...)))
end

struct RationalMap
	f::Function
	f_proj::Function

	df::Function
	df_proj::Function

	crit::Function

	function RationalMap(f, crit)
		if crit isa Function
            crit_function = to_point_map(crit)
        else
            crit_function = (_) -> Point(crit, 1.0)
        end

		df = derivative(f)
		return new(f, to_point_family(f), df, to_point_family(df), crit_function)
	end
end

function critical_points(func, parameter)
	@variables z, u, v, c

	f = func(z, c)
	df = expand_derivatives(Differential(z)(f))

	# Hack to write rational function as a ratio of polynomials
	# (Not sure it works on all cases.)
	df = df |>
		(x -> substitute(x, Dict(z => u / v))) |>
		simplify |>
		(x -> substitute(x, Dict(u => z, v => 1)))

	# Get numerator and denominator polynomials
	num, den = df |>
		Symbolics.value |>
		Symbolics.arguments

	num = substitute(num, Dict(c => parameter))
	den = substitute(den, Dict(c => parameter))

	num_coeffs = coeffs(num, z)
	den_coeffs = coeffs(den, z)

	num_roots = roots(Polynomial(num_coeffs))
	den_roots = roots(Polynomial(den_coeffs))

	pts = Point.(unique!(vcat(num_roots, den_roots)), 1)

	df_func = build_function(df, z, c, expression=Val{false})
	df_value_at_inf = to_point_family(df_func)(Point(1, 0), Point(parameter, 1))
	if isapprox(distance(df_value_at_inf, Point(1, 0)), 0) ||
			isapprox(distance(df_value_at_inf, Point(0, 1)), 0)
		push!(pts, Point(1, 0))
	end

	return pts
end

function coeffs(polynomial, z)
	higher_terms = polynomial
	coefficients = ComplexF64[]

	cutoff_test = 0
	while !isequal(higher_terms, 0) || cutoff_test > 10
		coefficient = substitute(higher_terms, Dict(z => 0))
		push!(coefficients, Symbolics.value(coefficient))

		higher_terms -= coefficient
		higher_terms = simplify(higher_terms / z)
		cutoff_test += 1
	end

	return coefficients
end