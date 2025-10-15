using Colors, ColorSchemes, FixedPointNumbers

const twilight_RGB = convert.(RGB24, ColorSchemes.twilight)
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

escape_time(f, z, c, a, ε, N) = to_color(convergence_time(f, z, c, a, ε, N)..., Val(:depth))
mod_period(f, z, c, a, ε, N) = to_color(convergence_time(f, z, c, a, ε, N)..., Val(:mod))
convergence_color(f, z, c, a, ε, N) = to_color(multiplier(f, z, c, a, ε, N))

struct ColoringData{T}
    method::Function
    attractors::Vector{Attractor{T}}
    update_attractors::Bool
end

get_attractor_type(::ColoringData{T}) where {T} = T
get_attractor_type(is_projective::Bool) = is_projective ? ProjectivePoint : ComplexF64

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
    return LCHab(50, chroma, hue), convert.(RGB24, gradient)
end

function generate_palette!(attractor::Attractor{T}, color) where {T<:PointLike}
    color = convert.(LCHab, color)
    chroma, hue = color.c, color.h

    palette = [LCHab(20 * cospi(t) + 50, chroma, hue) for t in range(0.0, 2.0, 510)]
    palette = convert.(RGB24, palette)

    attractor.color = color
    attractor.palette = palette

    return attractor
end

Attractor(c::T, m::Number, p::Integer) where {T} =
    Attractor{T}(c, 1, m, p, :twilight, twilight_RGB, true)
Attractor(c::Vector{T}, m::Number, p::Integer) where {T} =
    Attractor{T}(c, length(c), m, p, :twilight, twilight_RGB, true)

mutable struct ColoringScheme
    color::Union{Symbol,LCHab}
    palette::Vector{RGB24}
    continuous_coloring::Bool
end