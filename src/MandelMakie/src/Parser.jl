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
        df_dz = build_function(df_dz_expr, z, c, expression = Val{false})

        df_dc_expr = expand_derivatives(Differential(c)(expr))
        df_dc = build_function(df_dc_expr, z, c, expression = Val{false})

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
    return build_function(df, z, c, expression = Val{false})
end

function to_point_family(f)
    @variables u, v, c, d
    frac = simplify(f(u / v, c / d))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    fu = build_function(num, u, v, c, d, expression = Val{false})
    fv = build_function(den, u, v, c, d, expression = Val{false})

    return (pt_z, pt_c) -> normalize!(Point(fu(pt_z..., pt_c...), fv(pt_z..., pt_c...)))
end

function to_point_map(f)
    @variables c, d
    frac = simplify(f(c / d))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    fu = build_function(num, c, d, expression = Val{false})
    fv = build_function(den, c, d, expression = Val{false})

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

    # Hack to write rational function as a ratio of polynomials
    # (Not sure it works on all cases.)
    f =
        f |>
        (x -> substitute(x, Dict(z => u / v))) |>
        simplify |>
        (x -> substitute(x, Dict(u => z, v => 1)))

    # Get numerator and denominator polynomials
    p, q = f |> Symbolics.value |> Symbolics.arguments

    p = Polynomial(coeffs(substitute(p, Dict(c => parameter)), z))
    q = Polynomial(coeffs(substitute(q, Dict(c => parameter)), z))

    num = Polynomials.derivative(p) * q - p * Polynomials.derivative(q)
    num_roots = roots(num)

    pts = Point.(unique!(num_roots), 1)

    d = max(Polynomials.degree(p), Polynomials.degree(q))
    length(num_roots) < 2 * d - 2 && push!(pts, Point(1, 0))
    return pts
end

function coeffs(polynomial, z)
    higher_terms = polynomial
    coefficients = ComplexF64[]

    cutoff_test = 0
    while !isequal(higher_terms, 0) || cutoff_test > 20
        coefficient = substitute(higher_terms, Dict(z => 0))
        push!(coefficients, Symbolics.value(coefficient))

        higher_terms -= coefficient
        higher_terms = simplify(higher_terms / z)
        cutoff_test += 1
    end

    return coefficients
end
