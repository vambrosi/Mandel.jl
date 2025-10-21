using StaticArraysCore, LinearAlgebra, Symbolics

# --------------------------------------------------------------------------------------- #
# Complex Plane Definitions
# --------------------------------------------------------------------------------------- #

const ∞ = Inf + Inf * im

function distance(z::Number, w::Number)
    !isfinite(w) && return 1 / abs(z)
    return abs(z - w)
end

# --------------------------------------------------------------------------------------- #
# Complex Projective Line Definitions
# --------------------------------------------------------------------------------------- #

const ProjectivePoint = SVector{2,ComplexF64}

function Base.show(io::IO, pt::ProjectivePoint)
    z = pt[1] / pt[2]

    isfinite(z) ? print(io, z) : print(io, "∞")
end

Base.show(io::IO, ::MIME"text/plain", pt::ProjectivePoint) =
    print(io, "Point on the complex projective line:\n   ", pt)

Base.convert(::Type{ProjectivePoint}, z::ComplexF64) =
    isfinite(z) ? normalize(ProjectivePoint(z, 1)) : ProjectivePoint(1, 0)

function distance(pt1::ProjectivePoint, pt2::ProjectivePoint)
    return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

# --------------------------------------------------------------------------------------- #
# Change of Coordinates
# --------------------------------------------------------------------------------------- #

function projective_to_complex(pt::ProjectivePoint)
    z = pt[1] / pt[2]

    !isfinite(z) && return ∞
    return z
end

const PointLike = Union{ComplexF64,ProjectivePoint}
Base.convert(::Type{T}, z::Number) where {T<:PointLike} = ComplexF64(z)
Base.convert(::Type{ComplexF64}, z::ProjectivePoint) = projective_to_complex(z)

# --------------------------------------------------------------------------------------- #
# Extend functions to the Projective Line
# --------------------------------------------------------------------------------------- #

function extend_function(f::Function)
    # Makes sure function returns a vector
    if f(0.0im) isa Vector
        f_vector = z -> f(z)
    else
        f_vector = z -> [f(z)]
    end

    @variables u, v
    fractions = simplify.(f_vector(u / v))
    values = Symbolics.value.(fractions)

    # List of coordinate functions of the extended f from CP to CP x ... x CP
    functions = []

    for value in values
        if value isa Number
            # If expressions is a constant, push a constant function
            push!(functions, let
                pt = convert(ProjectivePoint, value)
                _ -> pt
            end)
        else
            num, den = Symbolics.arguments(value)
            fu = build_function(num, u, v, expression = Val{false})
            fv = build_function(den, u, v, expression = Val{false})

            push!(functions, z -> normalize(ProjectivePoint(fu(z[1], z[2]), fv(z[1], z[2]))))
        end
    end

    f_extended(z) = f_vector(z)
    f_extended(z::ProjectivePoint) = [g(z) for g in functions]
    return f_extended
end

function extend_family(f::Function)
    @variables u, v, a, b
    frac = simplify(f(u / v, a / b))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    f1 = build_function(num, u, v, a, b, expression = Val{false})
    f2 = build_function(den, u, v, a, b, expression = Val{false})

    g(z, c) = f(z, c)
    g(z::ProjectivePoint, c::Number) = normalize(ProjectivePoint(f1(z..., c, 1), f2(z..., c, 1)))
    g(z::ProjectivePoint, c::ProjectivePoint) = normalize(ProjectivePoint(f1(z..., c...), f2(z..., c...)))

    df_symb_00 = expand_derivatives(Differential(u)(simplify(f(u, a))))
    df_symb_0∞ = expand_derivatives(Differential(u)(simplify(1 / f(u, a))))
    df_symb_∞0 = expand_derivatives(Differential(u)(simplify(f(1 / u, a))))
    df_symb_∞∞ = expand_derivatives(Differential(u)(simplify(1 / f(1 / u, a))))

    df_00 = build_function(df_symb_00, u, a, expression = Val{false})
    df_0∞ = build_function(df_symb_0∞, u, a, expression = Val{false})
    df_∞0 = build_function(df_symb_∞0, u, a, expression = Val{false})
    df_∞∞ = build_function(df_symb_∞∞, u, a, expression = Val{false})

    function g(orbit::Vector{ComplexF64}, c::ComplexF64, ::Val{:diff})
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