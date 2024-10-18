module MandelGtk

export viewer

using Gtk4, Graphics, Cairo, StaticArraysCore, LinearAlgebra
# import Symbolics, Polynomials, Nemo
using Colors, ColorSchemes, FixedPointNumbers

function iterate(f, z, c, max_iter)
    for i in 1:max_iter
        z = f(z, c)

        if abs(z) > 10.0
            return i, abs(z)
        end
    end

    return nothing, abs(z)
end

to_color(::Nothing, ::Any) = get(ColorSchemes.twilight, 0.5)

function to_color(iterations::Integer, r::Real)
    height = (iterations - log(2, log(r))) / 50.0
    return get(ColorSchemes.twilight, mod(height, 1.0))
end

function to_complex_plane(roi, height, width, x, y)
    upp = roi.diameter / min(height, width)

    a = (x - width / 2) * upp
    b = -(y + 0.5 - height / 2) * upp

    return roi.center + complex(a, b)
end

@guarded function on_press(controller, n_press, x, y, view)
    view.pressed = true
    view.dragstart = [x, y]
    return
end

@guarded function on_release(controller, n_press, x, y, view, f)
    view.pressed || return
    view.pressed = false

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    z1 = to_complex_plane(view.roi, h, w, view.dragstart...)
    z2 = to_complex_plane(view.roi, h, w, x, y)

    view.roi.center += z1 - z2
    plt = plot(view.roi, f, w, h)

    surface = CairoImageSurface(plt)
    set_source_surface(ctx, surface)
    paint(ctx)
    reveal(canvas)

    view.plot = plt
    return
end

@guarded function on_motion(controller, x, y, view, coord)
    view.pointer = [x, y]

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    coord.label = string(to_complex_plane(view.roi, h, w, x, y))

    view.pressed || return

    set_source_rgb(ctx, 1, 1, 1)
    rectangle(ctx, 0, 0, w, h)
    fill(ctx)

    surface = CairoImageSurface(view.plot)

    Δx = x - view.dragstart[1]
    Δy = y - view.dragstart[2]

    set_source_surface(ctx, surface, Δx, Δy)
    paint(ctx)

    reveal(canvas)
    return
end

function zoom!(view, canvas, ctx, f, h, w, pointer)
    roi = deepcopy(view.roi)

    inv_scale = 1 / view.scale
    z = to_complex_plane(roi, h, w, view.pointer...)

    roi.diameter *= inv_scale
    roi.center = inv_scale * roi.center + (1 - inv_scale) * z

    plt = plot(roi, f, w, h)

    surface = CairoImageSurface(plt)
    set_source_surface(ctx, surface, 0, 0)
    paint(ctx)
    reveal(canvas)

    view.roi = roi
    view.scale = 1.0
    view.plot = plt

    return nothing
end

@guarded function on_zoom(controller, Δx, Δy, view, f)
    view.pressed && return
    view.zooming = true

    view.scale *= Δy < 0 ? 1.1 : 0.9
    view.scale = min(max(0.02, view.scale), 50)

    pointer = deepcopy(view.pointer)

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    set_source_rgb(ctx, 1, 1, 1)
    rectangle(ctx, 0, 0, w, h)
    fill(ctx)

    translate(ctx, pointer[1], pointer[2])
    scale(ctx, view.scale, view.scale)
    translate(ctx, -pointer[1], -pointer[2])
    surface = CairoImageSurface(view.plot)
    set_source_surface(ctx, surface, Δx, Δy)
    paint(ctx)
    reset_transform(ctx)
    reveal(canvas)

    close(view.zoom_timer)
    view.zoom_timer = Timer(0.1) do _
        zoom!(view, canvas, ctx, f, h, w, pointer)
    end

    return
end

mutable struct MandelROI
    center::ComplexF64
    diameter::Float64
end

function plot(mandel::MandelROI, f::Function, w::Integer, h::Integer)
    plt = Matrix{RGB24}(undef, w, h)
    futures = Vector{Task}(undef, h)

    upp = mandel.diameter / min(h, w)
    corner = mandel.center + upp * complex(-0.5 - w / 2, 1.0 + h / 2)

    for j in 1:h
        futures[j] = Threads.@spawn for i in 1:w
            c = corner + complex(i, -j) * upp
            plt[i, j] = to_color(iterate(f, 0.0im, c, 200)...)
        end
    end

    wait.(futures)
    return plt
end

mutable struct JuliaROI
    center::ComplexF64
    diameter::Float64
    parameter::ComplexF64
end

function plot(julia::JuliaROI, f::Function, w::Integer, h::Integer)
    plt = Matrix{RGB24}(undef, w, h)
    futures = Vector{Task}(undef, h)

    upp = julia.diameter / min(h, w)
    corner = julia.center + upp * complex(-0.5 - w / 2, 1.0 + h / 2)

    for j in 1:h
        futures[j] = Threads.@spawn for i in 1:w
            z = corner + complex(i, -j) * upp
            plt[i, j] = to_color(iterate(f, z, julia.parameter, 200)...)
        end
    end

    wait.(futures)
    return plt
end

const ROI = Union{MandelROI,JuliaROI}

mutable struct View
    roi::ROI
    plot::Matrix{RGB24}

    pressed::Bool
    zooming::Bool

    dragstart::Vector{Float64}
    pointer::Vector{Float64}
    scale::Float64

    zoom_timer::Timer
end

function View(roi::ROI, w::Integer, h::Integer)
    plt = Matrix{RGB24}(undef, w, h)

    pressed = false
    zooming = false

    dragstart = [0.0, 0.0]
    pointer = [0.0, 0.0]
    scale = 1.0

    timer = Timer(_ -> nothing, 0.1)
    return View(roi, plt, pressed, zooming, dragstart, pointer, scale, timer)
end

function add_mouse_events(canvas, view, f, coord_label)
    mouse_left = GtkGestureClick(canvas)
    signal_connect(mouse_left, "pressed") do args...
        on_press(args..., view)
    end

    signal_connect(mouse_left, "released") do args...
        on_release(args..., view, f)
    end

    mouse_motion = GtkEventControllerMotion(canvas)
    signal_connect(mouse_motion, "motion") do args...
        on_motion(args..., view, coord_label)
    end

    mouse_scroll =
        GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL, canvas)

    signal_connect(mouse_scroll, "scroll") do args...
        on_zoom(args..., view, f)
    end

    return mouse_left, mouse_motion, mouse_scroll
end

function add_mouse_events(mandel_canvas, julia_canvas, mandel, julia, f, coord_label)
    mouse_left = add_mouse_events(mandel_canvas, mandel, f, coord_label)[1]

    signal_connect(mouse_left, "released") do controller, n_press, x, y
        abs(mandel.dragstart[1] - x) > 5 && return
        abs(mandel.dragstart[2] - y) > 5 && return

        canvas = widget(controller)
        h = height(canvas)
        w = width(canvas)

        julia.roi.parameter = to_complex_plane(mandel.roi, h, w, x, y)

        h = height(julia_canvas)
        w = width(julia_canvas)

        ctx = getgc(julia_canvas)
        plt = plot(julia.roi, f, w, h)

        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)
        reveal(julia_canvas)

        julia.plot = plt
        return
    end
end

struct Viewer
    window::GtkWindow
    canvas::GtkCanvas
    overlay::GtkOverlay
end

function viewer(
    f;
    c = 0.0im,
    mandel_center = 0.0im,
    mandel_diameter = 4.0,
    julia_center = 0.0im,
    julia_diameter = 4.0,
    window_width = 800,
)
    win = GtkWindow("Mandel", window_width, window_width + 27, true, false)
    vbox = GtkBox(:v)
    push!(win, vbox)

    canvas = GtkCanvas()
    canvas.vexpand = canvas.hexpand = true

    overlay_size = window_width ÷ 4
    small_canvas = GtkCanvas(overlay_size, overlay_size)
    small_canvas.hexpand = small_canvas.vexpand = true

    frame = GtkFrame(small_canvas)
    frame.halign = 1
    frame.valign = 2
    frame.margin_start = 10
    frame.margin_bottom = 10

    overlay = GtkOverlay(canvas)
    add_overlay(overlay, frame)
    push!(vbox, overlay)

    coord_label = GtkLabel(string(0.0im))
    push!(vbox, coord_label)

    mandel = View(MandelROI(mandel_center, mandel_diameter), overlay_size, overlay_size)
    julia = View(JuliaROI(julia_center, julia_diameter, c), window_width, window_width)

    @guarded draw(canvas) do widget
        ctx = getgc(widget)
        h = height(widget)
        w = width(widget)

        d = min(h, w) ÷ 4
        small_canvas.content_height = d
        small_canvas.content_width = d

        plt = plot(julia.roi, f, w, h)
        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)

        julia.plot = plt
    end

    @guarded draw(small_canvas) do widget
        ctx = getgc(widget)
        h = height(widget)
        w = width(widget)

        plt = plot(mandel.roi, f, w, h)
        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)

        mandel.plot = plt
    end

    add_mouse_events(canvas, julia, f, coord_label)
    add_mouse_events(small_canvas, canvas, mandel, julia, f, coord_label)

    show(win)

    if !isinteractive()
        signal_connect(win, :close_request) do widget
            GC.gc()
        end
        @async Gtk4.GLib.glib_main()
    end

    return Viewer(win, canvas, overlay)
end

end