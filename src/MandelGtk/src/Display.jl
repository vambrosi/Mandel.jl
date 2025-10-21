using Gtk4, Gtk4.GLib, Graphics, Cairo

# --------------------------------------------------------------------------------------- #
# Main Definitions
# --------------------------------------------------------------------------------------- #

@enum PlaneType DynamicalPlane ParameterPlane

mutable struct ROI
    center::ComplexF64
    diameter::Float64
    parameter::ComplexF64
    plane_type::PlaneType
end

mutable struct Events
    pressed::Bool
    zooming::Bool

    dragstart::Vector{Float64}
    pointer::Vector{Float64}
    scale::Float64

    zoom_timer::Timer
end

mutable struct Options
    convergence_radius::Float64
    max_iterations::Int
    orbit_length::Int
    critical_length::Int
    compact_view::Bool
    is_family::Bool
    projective_metrics::Tuple{Bool,Bool}
    coloring_methods::Tuple{Symbol,Symbol}
    coloring_schemes::Vector{ColoringScheme}
end

mutable struct View
    roi::ROI
    init_roi::ROI

    canvas::GtkCanvas

    plot::Matrix{RGB24}
    coloring_data::ColoringData
    events::Events
end

struct Viewer
    window::GtkWindowLeaf
    mandel::View
    julia::View
end

# --------------------------------------------------------------------------------------- #
# Coordinate Change
# --------------------------------------------------------------------------------------- #

function canvas_to_complex(roi, height, width, x, y)
    upp = roi.diameter / min(height, width)

    a = (x - width / 2) * upp
    b = -(y - height / 2) * upp

    return roi.center + complex(a, b)
end

# --------------------------------------------------------------------------------------- #
# Options
# --------------------------------------------------------------------------------------- #

make_tuple(x) = (x, x)
make_tuple(x::Tuple) = x

function fix_criteria((mandel, julia))
    # Test if it is one of the valid options
    options = [:escape_time, :convergence_time, :mod_period, :preperiod]
    mandel in options || throw("Invalid Mandelbrot set coloring method")
    julia in options || throw("Invalid Julia set coloring method")

    # Mandelbrot can't use attractors so change to default
    mandel == :convergence_time && (mandel = :preperiod)
    mandel == :mod_period && (mandel = :preperiod)
    return mandel, julia
end

function get_coloring_data(map, c, coloring_method, projective_metric)
    if coloring_method == :escape_time
        method = escape_time
        a = get_attractor(map, c + 100, ∞, projective = projective_metric)
        attractors = [a]
        update_attractors = false
    elseif coloring_method == :convergence_time
        method = escape_time
        attractors = get_attractors(map, c, projective = projective_metric)
        update_attractors = true
    elseif coloring_method == :mod_period
        method = mod_period
        attractors = get_attractors(map, c, projective = projective_metric)
        update_attractors = true
    elseif coloring_method == :preperiod
        method = convergence_color
        attractors = Attractor{get_attractor_type(projective_metric)}[]
        update_attractors = false
    else
        throw("Invalid `coloring_method`")
    end

    return ColoringData(method, attractors, update_attractors)
end

function get_scheme(a::Attractor{T}) where {T}
    return ColoringScheme(a.color, a.palette, a.continuous_coloring)
end

function store_schemes!(options::Options, attractors::Vector{Attractor{T}}) where {T}
    options.coloring_schemes = [get_scheme(a) for a in attractors]
    return options
end

function sync_schemes!(options::Options, attractors::Vector{Attractor{T}}) where {T}
    schemes = options.coloring_schemes

    for (s, a) in zip(schemes, attractors)
        a.color = s.color
        a.palette = s.palette
        a.continuous_coloring = s.continuous_coloring
    end

    start = length(schemes) + 1

    for attractor in attractors[start:end]
        push!(schemes, get_scheme(attractor))
    end

    return options
end

# --------------------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------------------- #

function View(
    roi::ROI,
    canvas::GtkCanvas,
    coloring_data::ColoringData,
    w::Integer,
    h::Integer,
)
    plt = Matrix{RGB24}(undef, w, h)

    timer = Timer(_ -> nothing, 0.0)
    events = Events(false, false, [0.0, 0.0], [0.0, 0.0], 1.0, timer)

    return View(roi, deepcopy(roi), canvas, plt, coloring_data, events)
end

function parameter_slice!(plt, j, corner, w, upp, f, crit, coloring_data, options)
    attractor_type = get_attractor_type(coloring_data)

    for i in 1:w
        c = convert(attractor_type, corner + complex(i - 1, -j + 1) * upp)

        plt[i, j] = coloring_data.method(
            f,
            crit(c),
            c,
            coloring_data.attractors,
            options.convergence_radius,
            options.max_iterations,
        )
    end
end

function dynamical_slice!(plt, j, corner, w, upp, f, c, coloring_data, options)
    attractor_type = get_attractor_type(coloring_data)

    for i in 1:w
        z = convert(attractor_type, corner + complex(i - 1, -j + 1) * upp)

        plt[i, j] = coloring_data.method(
            f,
            z,
            c,
            coloring_data.attractors,
            options.convergence_radius,
            options.max_iterations,
        )
    end
end

function plot(
    roi::ROI,
    d_system::DynamicalSystem,
    options::Options,
    coloring_data::ColoringData,
    w::Integer,
    h::Integer,
)
    attractor_type = get_attractor_type(coloring_data)

    plt = Matrix{RGB24}(undef, w, h)
    futures = Vector{Task}(undef, h)

    upp = roi.diameter / min(h, w)
    corner = roi.center + upp * complex(0.5 - w / 2, h / 2 - 0.5)

    if roi.plane_type == ParameterPlane
        for j in 1:h
            futures[j] = Threads.@spawn parameter_slice!(
                plt,
                j,
                corner,
                w,
                upp,
                d_system.map,
                d_system.critical_point,
                coloring_data,
                options,
            )
        end
    else
        c = convert(attractor_type, roi.parameter)

        for j in 1:h
            futures[j] = Threads.@spawn dynamical_slice!(
                plt,
                j,
                corner,
                w,
                upp,
                d_system.map,
                c,
                coloring_data,
                options,
            )
        end
    end

    wait.(futures)
    return plt
end

function update_plot!(view, canvas, roi, d_system, options)
    Gtk4.GLib.g_idle_add() do
        w = width(canvas)
        h = height(canvas)
        ctx = getgc(canvas)

        plt = plot(roi, d_system, options, view.coloring_data, w, h)
        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)
        reveal(canvas)

        view.roi = roi
        view.events.scale = 1.0
        view.plot = plt

        false
    end
end

update_plot!(view, canvas, d_system, options) =
    update_plot!(view, canvas, view.roi, d_system, options)

# --------------------------------------------------------------------------------------- #
# Events
# --------------------------------------------------------------------------------------- #

@guarded function on_press(controller, n_press, x, y, view)
    view.events.pressed = true
    view.events.dragstart = [x, y]
    return
end

@guarded function on_release(controller, n_press, x, y, view, d_system, options)
    view.events.pressed || return
    view.events.pressed = false

    abs(view.events.dragstart[1] - x) < 10 &&
        abs(view.events.dragstart[2] - y) < 10 &&
        return

    canvas = widget(controller)

    h = height(canvas)
    w = width(canvas)

    z1 = canvas_to_complex(view.roi, h, w, view.events.dragstart...)
    z2 = canvas_to_complex(view.roi, h, w, x, y)

    view.roi.center += z1 - z2
    update_plot!(view, canvas, d_system, options)

    return
end

@guarded function on_motion(controller, x, y, view, coord)
    view.events.pointer = [x, y]

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    coord.label = string(canvas_to_complex(view.roi, h, w, x, y))

    view.events.pressed || return

    set_source_rgb(ctx, 1, 1, 1)
    rectangle(ctx, 0, 0, w, h)
    fill(ctx)

    surface = CairoImageSurface(view.plot)

    Δx = x - view.events.dragstart[1]
    Δy = y - view.events.dragstart[2]

    set_source_surface(ctx, surface, Δx, Δy)
    paint(ctx)

    reveal(canvas)
    return nothing
end

function zoom!(view, d_system, options, canvas, pointer)
    roi = deepcopy(view.roi)

    x, y = pointer
    inv_scale = 1 / view.events.scale

    h = height(canvas)
    w = width(canvas)

    # No idea why shifting by 1 is necessary here, but not elsewhere.
    z = canvas_to_complex(roi, h, w, x, y + 1)

    roi.diameter *= inv_scale
    roi.center = inv_scale * roi.center + (1 - inv_scale) * z

    update_plot!(view, canvas, roi, d_system, options)
    return nothing
end

@guarded function on_zoom(controller, Δx, Δy, view, d_system, options)
    view.events.pressed && return
    view.events.zooming = true

    view.events.scale *= Δy < 0 ? 1.1 : 0.9
    view.events.scale = min(max(0.02, view.events.scale), 50)

    pointer = deepcopy(view.events.pointer)

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    set_source_rgb(ctx, 1, 1, 1)
    rectangle(ctx, 0, 0, w, h)
    fill(ctx)

    translate(ctx, pointer[1], pointer[2])
    scale(ctx, view.events.scale, view.events.scale)
    translate(ctx, -pointer[1], -pointer[2])
    surface = CairoImageSurface(view.plot)
    set_source_surface(ctx, surface, Δx, Δy)
    paint(ctx)
    reset_transform(ctx)
    reveal(canvas)

    close(view.events.zoom_timer)
    view.events.zoom_timer = Timer(0.1) do _
        zoom!(view, d_system, options, canvas, pointer)
    end

    return
end

function add_events(canvas, view, d_system, options, coord_label)
    mouse_left = GtkGestureClick(canvas)
    signal_connect(mouse_left, "pressed") do args...
        on_press(args..., view)
    end

    signal_connect(mouse_left, "released") do args...
        on_release(args..., view, d_system, options)
    end

    mouse_motion = GtkEventControllerMotion(canvas)
    signal_connect(mouse_motion, "motion") do args...
        on_motion(args..., view, coord_label)
    end

    mouse_scroll =
        GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL, canvas)

    signal_connect(mouse_scroll, "scroll") do args...
        on_zoom(args..., view, d_system, options)
    end

    return mouse_left, mouse_motion, mouse_scroll
end

function add_events(
    mandel_canvas,
    julia_canvas,
    mandel,
    julia,
    d_system,
    options,
    coord_label,
)
    mouse_left = add_events(mandel_canvas, mandel, d_system, options, coord_label)[1]

    signal_connect(mouse_left, "released") do controller, n_press, x, y
        abs(mandel.events.dragstart[1] - x) > 10 && return
        abs(mandel.events.dragstart[2] - y) > 10 && return

        canvas = widget(controller)
        h = height(canvas)
        w = width(canvas)

        ctx = getgc(mandel_canvas)
        surface = CairoImageSurface(mandel.plot)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)
        reveal(mandel_canvas)

        julia.roi.parameter = canvas_to_complex(mandel.roi, h, w, x, y)

        coloring_data = julia.coloring_data
        if coloring_data.update_attractors
            T = get_attractor_type(coloring_data)
            attractor_list = get_attractors(
                d_system.map,
                julia.roi.parameter,
                projective = (T == ProjectivePoint),
                ε = options.convergence_radius / 1000,
            )
            julia.coloring_data = ColoringData{T}(
                coloring_data.method,
                attractor_list,
                coloring_data.update_attractors,
            )
            sync_schemes!(options, julia.coloring_data.attractors)
        end

        update_plot!(julia, julia_canvas, d_system, options)
        return
    end
end

# --------------------------------------------------------------------------------------- #
# Viewer
# --------------------------------------------------------------------------------------- #

function Viewer(
    f;
    crit = 0.0im,
    c = 0.0im,
    mandel_center = 0.0im,
    mandel_diameter = 4.0,
    julia_center = 0.0im,
    julia_diameter = 4.0,
    compact_view = true,
    coloring_method = :escape_time,
    projective_metric = false,
    window_width = 800,
)
    # Create Main Window
    win = GtkWindow("Mandel.jl", window_width, window_width + 15, true, false)
    vbox = GtkBox(:v)
    push!(win, vbox)

    julia_canvas = GtkCanvas()
    julia_canvas.vexpand = julia_canvas.hexpand = true

    overlay_size = window_width ÷ 4
    mandel_canvas = GtkCanvas(overlay_size, overlay_size)
    mandel_canvas.hexpand = mandel_canvas.vexpand = true

    is_mandel_small = true

    frame = GtkFrame(mandel_canvas)
    frame.halign = 1
    frame.valign = 2
    frame.margin_start = 10
    frame.margin_bottom = 10

    overlay = GtkOverlay(julia_canvas)
    add_overlay(overlay, frame)

    menu_button = GtkButton(; icon_name = "mail-send-receive-symbolic")
    menu_button.margin_start = 10
    menu_button.margin_top = 10
    menu_button.halign = 1
    menu_button.valign = 1

    add_overlay(overlay, menu_button)
    push!(vbox, overlay)

    coord_label = GtkLabel(string(0.0im))
    push!(vbox, coord_label)

    projective_metrics = make_tuple(projective_metric)
    coloring_methods = fix_criteria(make_tuple(coloring_method))

    is_family = hasmethod(f, Tuple{ComplexF64,ComplexF64})

    # Create Viewer Data
    d_system = DynamicalSystem(f, crit)
    options = Options(
        1e-3,
        300,
        1,
        1,
        compact_view,
        is_family,
        projective_metrics,
        coloring_methods,
        ColoringScheme[],
    )

    mandel_coloring =
        get_coloring_data(d_system.map, c, coloring_methods[1], projective_metrics[1])
    julia_coloring =
        get_coloring_data(d_system.map, c, coloring_methods[2], projective_metrics[2])

    mandel = View(
        ROI(mandel_center, mandel_diameter, 0.0im, ParameterPlane),
        mandel_canvas,
        mandel_coloring,
        overlay_size,
        overlay_size,
    )
    julia = View(
        ROI(julia_center, julia_diameter, c, DynamicalPlane),
        julia_canvas,
        julia_coloring,
        window_width,
        window_width,
    )

    # Add drawing handlers
    julia_timer = Timer(_ -> nothing, 0.0)
    @guarded draw(julia_canvas) do canvas
        w2 = width(canvas)
        h2 = height(canvas)
        ctx = getgc(canvas)

        set_source_rgb(ctx, 1, 1, 1)
        rectangle(ctx, 0, 0, w2, h2)
        fill(ctx)

        surface = CairoImageSurface(julia.plot)
        w1, h1 = size(julia.plot)

        set_source_surface(ctx, surface, (w2 - w1) ÷ 2, (h2 - h1) ÷ 2)
        paint(ctx)

        reveal(mandel_canvas)

        if is_mandel_small
            d = min(w2, h2) ÷ 4
            mandel_canvas.content_height = d
            mandel_canvas.content_width = d
        end

        close(julia_timer)
        julia_timer = Timer(0.1) do _
            update_plot!(julia, julia_canvas, d_system, options)
        end
    end

    mandel_timer = Timer(_ -> nothing, 0.0)
    @guarded draw(mandel_canvas) do canvas
        w2 = width(canvas)
        h2 = height(canvas)
        ctx = getgc(canvas)

        set_source_rgb(ctx, 1, 1, 1)
        rectangle(ctx, 0, 0, w2, h2)
        fill(ctx)

        surface = CairoImageSurface(mandel.plot)
        w1, h1 = size(mandel.plot)

        set_source_surface(ctx, surface, (w2 - w1) ÷ 2, (h2 - h1) ÷ 2)
        paint(ctx)

        reveal(mandel_canvas)

        if !is_mandel_small
            d = min(w2, h2) ÷ 4
            julia_canvas.content_height = d
            julia_canvas.content_width = d
        end

        close(mandel_timer)
        mandel_timer = Timer(0.1) do _
            update_plot!(mandel, mandel_canvas, d_system, options)
        end
    end

    # Add Events
    add_events(julia_canvas, julia, d_system, options, coord_label)
    add_events(mandel_canvas, julia_canvas, mandel, julia, d_system, options, coord_label)

    signal_connect(menu_button, "clicked") do widget
        frame_child = Gtk4.child(frame)
        overlay_child = Gtk4.child(overlay)

        h = height(overlay_child)
        w = width(overlay_child)
        d = min(h, w) ÷ 4

        overlay[] = nothing
        frame[] = nothing

        frame_child.content_height = 0
        frame_child.content_width = 0

        overlay_child.content_height = d
        overlay_child.content_width = d

        overlay[] = frame_child
        frame[] = overlay_child

        is_mandel_small = !is_mandel_small
    end

    show(win)

    if !isinteractive()
        signal_connect(win, :close_request) do widget
            GC.gc()
        end
        @async Gtk4.GLib.glib_main()
    end

    return Viewer(win, mandel, julia)
end

Base.show(io::IO, viewer::Viewer) = print(io, "MandelGtk Viewer")