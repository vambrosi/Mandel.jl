# --------------------------------------------------------------------------------------- #
# Parsing Functions
# --------------------------------------------------------------------------------------- #

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
	crit::Function

	function RationalMap(f, crit)
		if crit isa Function
            crit_function = to_point_map(crit)
        else
            crit_function = (_) -> Point(crit, 1.0)
        end

		return new(to_point_family(f), crit_function)
	end
end