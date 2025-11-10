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
    orbit_start::ComplexF64
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
    windows::Dict{Symbol,GtkWindowLeaf}
    mandel::View
    julia::View
    d_system::DynamicalSystem
    options::Options
end

const CLICK_DRIFT_TOLERANCE = 5

# --------------------------------------------------------------------------------------- #
# Coordinate Change
# --------------------------------------------------------------------------------------- #

function canvas_to_complex(roi, height, width, x, y)
    upp = roi.diameter / min(height, width)

    a = (x - width / 2) * upp
    b = -(y - height / 2) * upp

    return roi.center + complex(a, b)
end

function complex_to_canvas(roi, height, width, z::ComplexF64)
    ppu = min(height, width) / roi.diameter

    z -= roi.center

    x = width / 2 + ppu * real(z)
    y = height / 2 - ppu * imag(z)

    return x, y
end

function complex_to_canvas(roi, height, width, zs::Vector{ComplexF64})
    ppu = min(height, width) / roi.diameter

    zs .-= roi.center

    xs = real(zs) .* ppu .+ width / 2
    ys = imag(zs) .* (-ppu) .+ height / 2

    return collect(zip(xs, ys))
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

function redraw_orbit(ctx, roi, height, width, f, start, orbit_length, rgb)
    orbit = get_orbit(f, start, roi.parameter, orbit_length - 1)

    orbit_pxs = complex_to_canvas(roi, height, width, orbit)
    xy = orbit_pxs[1]

    set_source_rgb(ctx, rgb...)
    arc(ctx, xy..., 3, 0, 2pi)
    fill(ctx)

    set_line_width(ctx, 2)

    for iteration in 2:length(orbit_pxs)
        move_to(ctx, xy...)

        xy = orbit_pxs[iteration]
        line_to(ctx, xy...)
        stroke(ctx)

        arc(ctx, xy..., 3, 0, 2pi)
        fill(ctx)
    end

    return
end

function redraw_orbits(canvas, roi, d_system, options)
    w = width(canvas)
    h = height(canvas)
    ctx = getgc(canvas)

    if roi.plane_type == ParameterPlane
        set_source_rgb(ctx, 1, 0, 0)
        xy = complex_to_canvas(roi, h, w, roi.orbit_start)
        arc(ctx, xy..., 3, 0, 2pi)
        fill(ctx)
    else
        redraw_orbit(
            ctx,
            roi,
            h,
            w,
            d_system.map,
            roi.orbit_start,
            options.orbit_length,
            (1, 0, 0),
        )

        for c in d_system.critical_point(roi.parameter)
            redraw_orbit(
                ctx,
                roi,
                h,
                w,
                d_system.map,
                c,
                options.critical_length,
                (0, 0, 1),
            )
        end
    end
    return
end

function update_plot!(view, roi, d_system, options)
    @guarded Gtk4.GLib.g_idle_add() do
        w = width(view.canvas)
        h = height(view.canvas)
        ctx = getgc(view.canvas)

        plt = plot(roi, d_system, options, view.coloring_data, w, h)
        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)

        redraw_orbits(view.canvas, roi, d_system, options)

        reveal(view.canvas)

        view.roi = roi
        view.events.scale = 1.0
        view.plot = plt

        false
    end

    return
end

update_plot!(view, d_system, options) = update_plot!(view, view.roi, d_system, options)

# --------------------------------------------------------------------------------------- #
# Events
# --------------------------------------------------------------------------------------- #

@guarded function on_press(controller, n_press, x, y, view)
    view.events.pressed = true
    view.events.dragstart = [x, y]
    return
end

@guarded function on_release(controller, _, x, y, view, julia, d_system, options)
    view.events.pressed || return
    view.events.pressed = false

    canvas = widget(controller)

    h = height(canvas)
    w = width(canvas)

    if abs(view.events.dragstart[1] - x) < CLICK_DRIFT_TOLERANCE &&
       abs(view.events.dragstart[2] - y) < CLICK_DRIFT_TOLERANCE
        ctx = getgc(canvas)

        surface = CairoImageSurface(view.plot)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)

        view.roi.orbit_start = canvas_to_complex(view.roi, h, w, x, y)
        redraw_orbits(canvas, view.roi, d_system, options)
        reveal(canvas)

        if view.roi.plane_type == ParameterPlane
            julia.roi.parameter = canvas_to_complex(view.roi, h, w, x, y)

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

            update_plot!(julia, d_system, options)
        end
        return
    end

    z1 = canvas_to_complex(view.roi, h, w, view.events.dragstart...)
    z2 = canvas_to_complex(view.roi, h, w, x, y)

    view.roi.center += z1 - z2
    update_plot!(view, d_system, options)

    return
end

@guarded function on_motion(controller, x, y, view, d_system, options, coord)
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

    translate(ctx, Δx, Δy)
    redraw_orbits(canvas, view.roi, d_system, options)
    translate(ctx, -Δx, -Δy)

    reveal(canvas)
    return
end

function zoom!(view, d_system, options, canvas, pointer)
    roi = deepcopy(view.roi)

    x, y = pointer
    inv_scale = 1 / view.events.scale

    h = height(canvas)
    w = width(canvas)

    z = canvas_to_complex(roi, h, w, x, y)

    roi.diameter *= inv_scale
    roi.center = inv_scale * roi.center + (1 - inv_scale) * z

    update_plot!(view, roi, d_system, options)
    return
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
    set_source_surface(ctx, surface, 0, 0)
    paint(ctx)

    redraw_orbits(canvas, view.roi, d_system, options)
    reset_transform(ctx)
    reveal(canvas)

    close(view.events.zoom_timer)
    view.events.zoom_timer = Timer(0.1) do _
        zoom!(view, d_system, options, canvas, pointer)
    end

    return
end

function add_events(canvases, views, d_system, options, coord_label)
    julia = views[2].roi.plane_type == DynamicalPlane ? views[2] : views[1]

    for (canvas, view) in zip(canvases, views)
        mouse_left = GtkGestureClick(canvas)
        signal_connect(mouse_left, "pressed") do args...
            on_press(args..., view)
        end

        signal_connect(mouse_left, "released") do args...
            on_release(args..., view, julia, d_system, options)
        end

        mouse_motion = GtkEventControllerMotion(canvas)
        signal_connect(mouse_motion, "motion") do args...
            on_motion(args..., view, d_system, options, coord_label)
        end

        mouse_scroll =
            GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL, canvas)

        signal_connect(mouse_scroll, "scroll") do args...
            on_zoom(args..., view, d_system, options)
        end
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
    c = ComplexF64(c)

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

    button_box = GtkBox(:h)
    button_box.margin_start = 10
    button_box.margin_top = 10
    button_box.halign = 1
    button_box.valign = 1

    menu = GMenu()
    item2 = GMenuItem("Options...", "win.options")
    push!(menu, item2)
    item3 = GMenuItem("Quit", "win.close")
    push!(menu, item3)

    menu_button = GtkMenuButton(; icon_name = "open-menu-symbolic")
    Gtk4.menu_model(menu_button, menu)
    menu_button.tooltip_text = "Menu"
    menu_button.halign = 0
    menu_button.valign = 1

    switch_button = GtkButton(; icon_name = "mail-send-receive-symbolic")
    switch_button.tooltip_text = "Switch Views"

    reset_button = GtkButton(; icon_name = "zoom-fit-best-symbolic")
    reset_button.tooltip_text = "Reset View"
    push!(button_box, menu_button, switch_button, reset_button)

    add_overlay(overlay, button_box)
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
        ROI(mandel_center, mandel_diameter, 0.0im, ParameterPlane, c),
        mandel_canvas,
        mandel_coloring,
        overlay_size,
        overlay_size,
    )

    julia = View(
        ROI(julia_center, julia_diameter, c, DynamicalPlane, julia_center),
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
            update_plot!(julia, d_system, options)
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
            update_plot!(mandel, d_system, options)
        end
    end

    # Add Events
    add_events(
        [mandel_canvas, julia_canvas],
        [mandel, julia],
        d_system,
        options,
        coord_label,
    )

    @guarded signal_connect(switch_button, "clicked") do widget
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
        return
    end

    signal_connect(reset_button, "clicked") do widget
        overlay_child = Gtk4.child(overlay)

        view = overlay_child == julia_canvas ? julia : mandel
        update_plot!(view, view.init_roi, d_system, options)
        return
    end

    windows = Dict(:main => win)

    # Menu actions
    function options_action(a, par)::Nothing
        if haskey(windows, :options)
            present(windows[:options])
            return
        end

        dialog = GtkWindow("Options", 200, 100, true, false)
        windows[:options] = dialog

        content = GtkBox(:v)
        content.margin_start = 10
        content.margin_end = 10
        content.margin_top = 10
        content.margin_bottom = 10
        content.spacing = 10

        max_label = GtkLabel("Max Iterations:")
        max_spin = GtkSpinButton(1, 10000, 1)
        max_spin.value = options.max_iterations

        point_label = GtkLabel("Point Orbit Length:")
        point_spin = GtkSpinButton(1, 10000, 1)
        point_spin.value = options.orbit_length

        critical_label = GtkLabel("Critical Orbit Length:")
        critical_spin = GtkSpinButton(1, 10000, 1)
        critical_spin.value = options.critical_length

        convergence_label = GtkLabel("Convergence Radius:")
        convergence_spin = GtkEntry()
        convergence_spin.text = string(options.convergence_radius)

        button = GtkButton(; label = "Update Values")

        push!(content, max_label)
        push!(content, max_spin)
        push!(content, point_label)
        push!(content, point_spin)
        push!(content, critical_label)
        push!(content, critical_spin)
        push!(content, convergence_label)
        push!(content, convergence_spin)
        push!(content, button)
        push!(dialog, content)

        show(dialog)

        signal_connect(button, "clicked") do widget
            options.max_iterations = Int(max_spin.value)
            options.orbit_length = Int(point_spin.value)
            options.critical_length = Int(critical_spin.value)
            options.convergence_radius = parse(Float64, convergence_spin.text)

            update_plot!(mandel, d_system, options)
            update_plot!(julia, d_system, options)
        end

        signal_connect(dialog, :close_request) do widget
            delete!(windows, :options)
            return
        end

        return
    end

    function close_action(a, par)::Nothing
        haskey(windows, :options) && destroy(windows[:options])
        destroy(windows[:main])
        return
    end

    action_group = GSimpleActionGroup()
    add_action(GActionMap(action_group), "close", close_action)
    add_action(GActionMap(action_group), "options", options_action)
    push!(win, Gtk4.GLib.GActionGroup(action_group), "win")

    signal_connect(win, :close_request) do widget
        haskey(windows, :options) && destroy(windows[:options])
        return
    end

    show(win)

    if !isinteractive()
        @async Gtk4.GLib.glib_main()
        Gtk4.GLib.waitforsignal(win, :close_request)
    end

    return Viewer(windows, mandel, julia, d_system, options)
end

Base.show(io::IO, viewer::Viewer) = print(io, "MandelGtk Viewer")