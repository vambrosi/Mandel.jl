module MandelGtk

export viewer

using Gtk4, Graphics, Cairo, Colors, FixedPointNumbers

include("ColorSchemes.jl")

function iterate(f, z, c, max_iter)
    for i in 1:max_iter
        z = f(z, c)

        if abs(z) > 10.0
            return i, abs(z)
        end
    end

    return nothing, abs(z)
end

color(::Nothing, ::Any) = twilight_RGB24[500]

function color(iterations::Integer, r::Real)
    height = (iterations - log(2, log(r))) / 50.0
    index = mod1(round(Int, 999 * mod(height, 1.0)), 1000)
    return twilight_RGB24[index]
end

function to_complex_plane(box, height, width, x, y)
    upp = box.diameter / min(height, width)

    a = (x + 1 - width / 2) * upp
    b = -(y + 2 - height / 2) * upp

    return box.center + complex(a, b)
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

    z1 = to_complex_plane(view.box, h, w, view.dragstart...)
    z2 = to_complex_plane(view.box, h, w, x, y)

    view.box.center += z1 - z2
    plt = plot(view.box, f, w, h)

    surface = CairoImageSurface(plt)
    set_source_surface(ctx, surface)
    paint(ctx)
    reveal(canvas)

    view.plot = plt
    return
end

@guarded function on_motion(controller, x, y, view)
    view.pointer = [x, y]
    view.pressed || return

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

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
    box = deepcopy(view.box)

    inv_scale = 1 / view.scale
    z = to_complex_plane(box, h, w, view.pointer...)

    box.diameter *= inv_scale
    box.center = inv_scale * box.center + (1 - inv_scale) * z

    plt = plot(box, f, w, h)

    surface = CairoImageSurface(plt)
    set_source_surface(ctx, surface, 0, 0)
    paint(ctx)
    reveal(canvas)

    view.box = box
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

mutable struct MandelBox
    center::ComplexF64
    diameter::Float64
end

mutable struct JuliaBox
    center::ComplexF64
    diameter::Float64
    parameter::ComplexF64
end

mutable struct View
    box::Union{MandelBox,JuliaBox}
    plot::Matrix{RGB24}

    pressed::Bool
    zooming::Bool

    dragstart::Vector{Float64}
    pointer::Vector{Float64}
    scale::Float64

    zoom_timer::Timer
end

function plot(mandel::MandelBox, f::Function, w::Integer, h::Integer)
    plt = Matrix{RGB24}(undef, w, h)
    futures = Vector{Task}(undef, h)

    upp = mandel.diameter / min(h, w)
    corner = mandel.center + upp * complex(0.5 - w / 2, -0.5 + h / 2)

    for j in 1:h
        futures[j] = Threads.@spawn for i in 1:w
            c = corner + complex(i, -j) * upp
            plt[i, j] = color(iterate(f, 0.0im, c, 200)...)
        end
    end

    wait.(futures)
    return plt
end

function plot(julia::JuliaBox, f::Function, w::Integer, h::Integer)
    plt = Matrix{RGB24}(undef, w, h)
    futures = Vector{Task}(undef, h)

    upp = julia.diameter / min(h, w)
    corner = julia.center + upp * complex(0.5 - w / 2, -0.5 + h / 2)

    for j in 1:h
        futures[j] = Threads.@spawn for i in 1:w
            z = corner + complex(i, -j) * upp
            plt[i, j] = color(iterate(f, z, julia.parameter, 200)...)
        end
    end

    wait.(futures)
    return plt
end

struct Viewer
    window::GtkWindow
    canvas::GtkCanvas
end

function viewer(f::Function, c::Number, is_mandel::Bool)
    win = GtkWindow("Mandel", 500, 500, true, false)
    vbox = GtkBox(:v)
    push!(win, vbox)

    canvas = GtkCanvas()
    canvas.vexpand = canvas.hexpand = true
    push!(vbox, canvas)

    plt = Matrix{RGB24}(undef, 500, 500)

    pressed = false
    zooming = false

    dragstart = [0.0, 0.0]
    pointer = [0.0, 0.0]
    scale = 1.0

    box = is_mandel ? MandelBox(-0.5, 4.0) : JuliaBox(0.0im, 4.0, c)
    view = View(
        box,
        plt,
        pressed,
        zooming,
        dragstart,
        pointer,
        scale,
        Timer(_ -> nothing, 0.1),
    )

    @guarded draw(canvas) do widget
        ctx = getgc(canvas)
        h = height(canvas)
        w = width(canvas)

        plt = plot(view.box, f, w, h)
        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)

        view.plot = plt
    end

    mouse_left = GtkGestureClick(canvas)
    signal_connect(mouse_left, "pressed") do args...
        on_press(args..., view)
    end

    signal_connect(mouse_left, "released") do args...
        on_release(args..., view, f)
    end

    mouse_motion = GtkEventControllerMotion(canvas)
    signal_connect(mouse_motion, "motion") do args...
        on_motion(args..., view)
    end

    mouse_scroll =
        GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL, canvas)

    signal_connect(mouse_scroll, "scroll") do args...
        on_zoom(args..., view, f)
    end

    show(win)

    if !isinteractive()
        signal_connect(win, :close_request) do widget
            GC.gc()
        end
        @async Gtk4.GLib.glib_main()
    end

    return Viewer(win, canvas)
end

end