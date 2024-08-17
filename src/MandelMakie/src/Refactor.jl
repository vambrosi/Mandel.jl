module Refactor

export Viewer

using GLMakie, Symbolics, StaticArraysCore, LinearAlgebra, Polynomials
import Dates

# --------------------------------------------------------------------------------------- #
# Complex Plane Definitions
# --------------------------------------------------------------------------------------- #

const ∞ = NaN + NaN * im

function distance(z::Number, w::Number)
    !isfinite(w) && return 1 / abs(z)
    return abs(z - w)
end

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

const Point = SVector{2,ComplexF64}

to_projective(z) = normalize(Point(z, 1))

function distance(pt1::Point, pt2::Point)
    return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

function extend_function(f::Function)
    @variables u, v
    frac = simplify(f(u / v))
    value = Symbolics.value(frac)

    if value isa Number
        h = let pt = to_projective(value)
            _ -> pt
        end
    else
        num, den = Symbolics.arguments(value)
        fu = build_function(num, u, v, expression = Val{false})
        fv = build_function(den, u, v, expression = Val{false})

        h = z -> normalize(Point(fu(z...), fv(z...)))
    end

    g(z::ComplexF64) = f(z)
    g(z::Point) = h(z)
    return g
end

function extend_family(f::Function)
    @variables u, v, c, d
    frac = simplify(f(u / v, c / d))
    value = Symbolics.value(frac)

    num, den = Symbolics.arguments(value)
    fu = build_function(num, u, v, c, d, expression = Val{false})
    fv = build_function(den, u, v, c, d, expression = Val{false})

    g(z::Number, c::Number) = f(z, c)
    g(z::Point, c::Number) = normalize(Point(fu(z..., c, 1), fv(z..., c, 1)))
    g(z::Point, c::Point) = normalize(Point(fu(z..., c...), fv(z..., c...)))

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

        return new(extend_family(h), extend_function(critical_point))
    end
end

function DynamicalSystem(f::Function, c::Number)
    c = convert(ComplexF64, c)
    g = let c = c
        _ -> c
    end

    return DynamicalSystem(f, g)
end

struct Attractor
    cycle::Union{ComplexF64,Vector{ComplexF64}}
    multiplier::ComplexF64
    power::Int
end

const MaybeAttractor = Union{Nothing,Attractor}

const empty_attractor = nothing

function iscloseby(z, w::ComplexF64, ε)
    distance(z, w) < ε && return true

    return false
end

function iscloseby(z, cycle::Vector, ε)
    for w in cycle
        distance(z, w) < ε && return true
    end

    return false
end

function convergence_time(f, z, c, attractors, ε, max_iterations)
    for iteration = 0:max_iterations
        for attractor in attractors
            if iscloseby(z, attractor.cycle, ε)
                return iteration, z, attractor
            end
        end

        z = f(z, c)
    end

    return max_iterations + 1, z, empty_attractor
end

const PointLike = Union{ComplexF64,Point}
Base.convert(::Type{PointLike}, z::Number) = ComplexF64(z)

struct CloseBy
    iterations::Int
    z::PointLike
    w::PointLike
end

const MaybeCloseBy = Union{Nothing,CloseBy}

function period_multiple_apart(f, z, c, ε, max_iterations)::MaybeCloseBy
    slow = fast = z

    for n = 1:max_iterations
        slow = f(slow, c)
        fast = f(f(fast, c), c)

        distance(slow, fast) <= ε && return CloseBy(n, slow, fast)
    end

    return nothing
end

function iterate_until_close(f, z, w, c, ε, max_iterations)::MaybeCloseBy
    for n = 1:max_iterations
        z = f(z, c)

        distance(z, w) < ε && return CloseBy(n, z, w)
    end

    return nothing
end

function iterate_both_until_close(f, z, w, c, ε, max_iterations)::MaybeCloseBy
    for n = 1:max_iterations
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

function multiplier(f, z0, c, _, ε, max_iterations)::MaybeOrbitData
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

# --------------------------------------------------------------------------------------- #
# Coloring Algorithms
# --------------------------------------------------------------------------------------- #

const DEFAULT_COLOR_LEVEL = 0.5

function to_color(approach::OrbitData)
    return mod(approach.preperiod / approach.period / 64.0, 1.0)
end

function to_color(iterations::Integer, z::Number, attractor::MaybeAttractor)
    attractor == empty_attractor && return DEFAULT_COLOR_LEVEL
    mod((iterations + 1.0 - log(attractor.power, log(abs(z)))) / 64.0, 1.0)
end

to_color(::Nothing) = DEFAULT_COLOR_LEVEL

escape_time(f, z, c, a, ε, N) = to_color(convergence_time(f, z, c, a, ε, N)...)
convergence_color(f, z, c, a, ε, N) = to_color(multiplier(f, z, c, a, ε, N))

struct ColoringData
    method::Function
    is_projective::Bool
    attractors::Vector{Attractor}
end

# --------------------------------------------------------------------------------------- #
# Views of Mandelbrot and Julia Sets
# --------------------------------------------------------------------------------------- #

abstract type View end

mutable struct Options
    convergence_radius::Float64
    max_iterations::Int
    orbit_length::Int
    critical_length::Int
    compact_view::Bool
    is_family::Bool
end

mutable struct MandelView <: View
    center::ComplexF64
    diameter::Float64
    pixels::Int

    init_center::ComplexF64
    init_diameter::Float64

    color_levels::Observable{Matrix{Float64}}
    points::Observable{Vector{ComplexF64}}
    marks::Observable{Vector{ComplexF64}}

    coloring_data::ColoringData

    function MandelView(; center, diameter, pixels, coloring_data)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        color_levels = zeros(Float64, pixels, pixels)
        points = ComplexF64[center]
        marks = ComplexF64[]

        return new(
            center,
            diameter,
            pixels,
            center,
            diameter,
            color_levels,
            points,
            marks,
            coloring_data,
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

    color_levels::Observable{Matrix{Float64}}
    points::Observable{Vector{ComplexF64}}
    marks::Observable{Vector{ComplexF64}}

    coloring_data::ColoringData

    function JuliaView(; center, diameter, parameter, pixels, coloring_data)
        diameter > 0.0 || throw("diameter must be a positive real")
        pixels > 0 || throw("pixels must be a positive integer")

        color_levels = zeros(Float64, pixels, pixels)
        points = ComplexF64[center]
        marks = ComplexF64[]

        return new(
            center,
            diameter,
            parameter,
            pixels,
            center,
            diameter,
            color_levels,
            points,
            marks,
            coloring_data,
        )
    end
end

# --------------------------------------------------------------------------------------- #
# Change of Coordinates between Complex Plane and Pixel Space
# --------------------------------------------------------------------------------------- #

function to_complex_plane(view::View, pixel_vector)
    x, y = pixel_vector[1], pixel_vector[2]

    upp = view.diameter / view.pixels
    a = (x + 0.5 - view.pixels / 2) * upp
    b = (y + 0.5 - view.pixels / 2) * upp

    return view.center + complex(a, b)
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
    @inbounds for i = 1:pxs
        c =
            coloring_data.is_projective ? to_projective(corner + step * complex(i, j)) :
            corner + step * complex(i, j)

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

    @inbounds for j = 1:view.pixels
        futures[j] = Threads.@spawn mandel_slice!(
            view.color_levels[],
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
    @inbounds for i = 1:pxs
        z =
            coloring_data.is_projective ? to_projective(corner + step * complex(i, j)) :
            corner + step * complex(i, j)

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
    parameter =
        view.coloring_data.is_projective ? to_projective(view.parameter) : view.parameter
    futures = Vector{Task}(undef, view.pixels)

    @inbounds for j = 1:view.pixels
        futures[j] = Threads.@spawn julia_slice!(
            view.color_levels[],
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
    notify(view.color_levels)
    notify(view.points)
    notify(view.marks)

    return view
end

function pick_parameter!(
    julia::JuliaView,
    mandel::MandelView,
    d_system::DynamicalSystem,
    options::Options,
    point,
)
    julia.parameter = point

    julia.marks[] = orbit(
        d_system.map,
        d_system.critical_point(julia.parameter),
        julia.parameter,
        options.critical_length - 1,
    )

    mandel.points[] = [julia.parameter]

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
    end

    for axis in [left_axis, right_axis]
        hidedecorations!(axis)
        deregister_interaction!(axis, :rectanglezoom)
        deregister_interaction!(axis, :limitreset)
        deregister_interaction!(axis, :dragpan)
    end

    return Frame(left_axis, mandel, Dict()), Frame(right_axis, julia, Dict())
end

function create_plot!(frame::Frame)
    empty!(frame.axis)
    view = frame.view[]

    heatmap!(
        frame.axis,
        view.color_levels,
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

    mark_vectors = lift(view.marks) do zs
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

    heatmap!(ax, view.color_levels[], colormap = :twilight, colorrange = (0.0, 1.0))

    Makie.save(filename, fig)
end

function add_frame_events!(
    frame::Frame,
    topframe::Frame,
    d_system::DynamicalSystem,
    options::Options,
    julia::JuliaView,
)
    axis = frame.axis
    scene = axis.scene

    # Mouse Events
    dragging = false
    dragstart = Point2f(0.0)
    dragend = Point2f(0.0)

    to_world_at_start = z -> to_world(scene, z)

    is_zooming = false
    zooming = Timer(_ -> zoom!(frame, d_system, options), 0.1)

    is_topframe = frame == topframe
    z_level = is_topframe ? 10 : 0

    on(events(scene).mousebutton) do event
        is_zooming && return Consume(false)

        mp = mouseposition_px(scene)
        view = frame.view[]

        if event.button == Mouse.right
            if event.action == Mouse.press &&
               is_mouseinside(axis) &&
               (is_topframe || !is_mouseinside(topframe.axis))
                dragging = true
                to_world_at_start = z -> to_world(scene, z)
                dragstart = to_world_at_start(mp)

            elseif event.action == Mouse.release && dragging
                dragend = to_world_at_start(mp)
                view.center +=
                    to_complex_plane(view, dragstart) - to_complex_plane(view, dragend)
                translate!(scene, 0, 0, z_level)
                update_view!(view, d_system, options)
                reset_limits!(axis)
                dragging = false
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
                end

                point = to_complex_plane(view, to_world_at_start(mp))
                if view isa MandelView
                    pick_parameter!(julia, view, d_system, options, point)
                elseif view isa JuliaView
                    pick_orbit!(view, d_system, options, point)
                end
            end
        end

        return Consume(false)
    end

    on(events(scene).mouseposition) do event
        if dragging
            mp = mouseposition(scene)
            translate!(scene, mp - dragstart..., z_level)
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

    frame.events[:dragging] = dragging
    frame.events[:dragstart] = dragstart
    frame.events[:dragend] = dragstart
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
    end
end

function add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)
    layout = GridLayout(figure[2, 1], tellwidth = false)
    button_shift = options.is_family ? 2 : 0

    labels = Dict(
        :max_iter => Label(layout[1, button_shift+1], "Maximum\nIterations:"),
        :orbit_len => Label(layout[1, button_shift+3], "Orbit\nLength:"),
        :critical_length =>
            Label(layout[1, button_shift+5], "Critical Point\nOrbit Length:"),
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
        julia.marks[] = orbit(
            d_system.map,
            d_system.critical_point(julia.parameter),
            julia.parameter,
            options.critical_length - 1,
        )
    end

    on(inputs[:convergence_radius].stored_string) do s
        options.convergence_radius = parse(Float64, s)
        update_view!(julia, d_system, options)
        options.is_family && update_view!(mandel, d_system, options)
    end

    if options.is_family
        inputs[:switch_layout] = Button(layout[1, 1], label = "↰", halign = :left)
        inputs[:switch_positions] = Button(layout[1, 2], label = "↔", halign = :left)

        on(inputs[:switch_layout].clicks, priority = 300) do event
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

        on(inputs[:switch_positions].clicks, priority = 300) do event
            left_frame.view[], right_frame.view[] = right_frame.view[], left_frame.view[]
            create_plot!(left_frame)
            create_plot!(right_frame)
        end
    end

    return inputs
end

"""
	Viewer(f; <keyword arguments>)

Create a `Viewer` that plots the Mandelbrot and the Julia sets associated with the \
function `f`. `f` can have one or two inputs, the second being the parameter. \
By default, the `Viewer` opens as a separate window.

# Examples
```julia-repl
julia> Viewer((z, c) -> z^2 + c, mandel_center=-0.5)
```
```julia-repl
julia> f(z, λ) = z^2 + λ / z^2
julia> crit(λ) = λ^(1/4)
julia> Viewer(f; crit=crit, mandel_diameter=1.0)
```

# Shortcuts

- `Right Click Drag`: Pans view.
- `Mouse Scroll`: Zooms in or out.
- `Left Click` (parameter space): Chooses parameter.
- `Left Click` (dynamical space): Chooses orbit initial point.
- `Ctrl + Left Click`: Resets view to initial `center` and `diameter`.
- `Ctrl + S`: Saves view that is under the mouse pointer (files in `./imgs`).

# Arguments
- `crit = 0.0im`: Function that gives a critical point for each parameter. Used to \
to plot the Mandelbrot set. If it is a constant function, you can just input the constant \
directly.

- `c = 0.0im`: Initial parameter used to plot the Julia set.

- `mandel_center = 0.0im`: Initial center of the Mandelbrot plot.

- `mandel_diameter = 4.0`: Initial diameter of the Mandelbrot plot.

- `julia_center = 0.0im`: Initial center of the Julia plot.

- `julia_diameter = 4.0`: Initial diameter of the Julia plot.

- `grid_width = 800`: Width (and height) of the grid of complex numbers used to plot sets.

- `compact_view = true`: If 'true' one of the plots is show as an inset plot, if `false` \
they are shown side-by-side.

- `mandel_coloring_method = :escape_time`: Chooses the coloring method for the Mandelbrot \
plot. The options are `:escape_time`, `:plane_convergence`, `:projective_convergence`.

- `julia_coloring_method = :escape_time`: Same as for `mandel` version, but with two \
extra options: `:plane_limiting_attractor`, and `:projective_limiting_attractor`.
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
        mandel_coloring_method = :escape_time,
        julia_coloring_method = :escape_time,
    )

        is_family = hasmethod(f, Tuple{ComplexF64,ComplexF64})

        coloring_data = Dict(
            :escape_time => ColoringData(escape_time, false, [Attractor(∞, 0, 2)]),
            :plane_convergence => ColoringData(convergence_color, false, Attractor[]),
            :projective_convergence =>
                ColoringData(convergence_color, false, Attractor[]),
        )

        d_system = DynamicalSystem(f, crit)
        options = Options(1e-3, 200, 1, 1, compact_view, is_family)
        figure = Figure(size = (750, 650))

        mandel =
            !is_family ? nothing :
            MandelView(
                center = mandel_center,
                diameter = mandel_diameter,
                pixels = grid_width,
                coloring_data = coloring_data[mandel_coloring_method],
            )

        julia = JuliaView(
            center = julia_center,
            diameter = julia_diameter,
            parameter = c,
            pixels = grid_width,
            coloring_data = coloring_data[julia_coloring_method],
        )

        julia.marks[] = [d_system.critical_point(julia.parameter)]
        pick_orbit!(julia, d_system, options, julia.points[][begin])

        left_frame, right_frame = create_frames!(figure, options, mandel, julia)

        if is_family
            add_frame_events!(right_frame, left_frame, d_system, options, julia)
            add_frame_events!(left_frame, left_frame, d_system, options, julia)

            create_plot!(left_frame)
            update_view!(mandel, d_system, options)
            mandel.points[] = [julia.parameter]
        else
            add_frame_events!(right_frame, right_frame, d_system, options, julia)
        end

        create_plot!(right_frame)
        update_view!(julia, d_system, options)

        colgap!(content(figure[1, 1]), 10)
        inputs =
            add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)

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

end
