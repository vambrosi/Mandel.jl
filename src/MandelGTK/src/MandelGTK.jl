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

function to_complex_plane(view, height, width, x, y)
    upp = view.diameter / min(height, width)

    a = (x + 0.5 - width / 2) * upp
    b = -(y + 0.5 - height / 2) * upp

    return view.center + complex(a, b)
end

function on_release(controller, n_press, x, y, view, f)
    view.pressed || return
    view.pressed = false

    canvas = widget(controller)
    ctx = getgc(canvas)

    h = height(canvas)
    w = width(canvas)

    Gtk4.GLib.g_idle_add() do
        z1 = to_complex_plane(view, h, w, view.dragstart...)
        z2 = to_complex_plane(view, h, w, x, y)

        view.center += z1 - z2
        plt = plot!(view, f, w, h)

        surface = CairoImageSurface(plt)
        set_source_surface(ctx, surface)
        paint(ctx)
        reveal(canvas)
        false
    end

    return
end

function on_motion(controller, x, y, view)
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

# function on_zoom(controller, x, y, view, f)
#     view.pressed && return

#     canvas = widget(controller)
#     ctx = getgc(canvas)

#     h = height(canvas)
#     w = width(canvas)

#     set_source_rgb(ctx, 1, 1, 1)
#     rectangle(ctx, 0, 0, w, h)
#     fill(ctx)

#     translate(ctx, view.pointer[1], view.pointer[2])
#     scale(ctx, view.scale, view.scale)
#     paint(ctx)
#     translate(ctx, 0, 0)

#     reveal(canvas)
# end

mutable struct MandelView
    center::ComplexF64
    diameter::Float64
    plot::Matrix{RGB24}

    pressed::Bool
    zooming::Bool

    dragstart::Vector{Float64}
    pointer::Vector{Float64}
    scale::Float64
end

function plot!(mandel::MandelView, f::Function, w::Integer, h::Integer)
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
    mandel.plot = plt

    return plt
end

mutable struct JuliaView
    center::ComplexF64
    diameter::Float64
    parameter::ComplexF64
    plot::Matrix{RGB24}

    pressed::Bool
    zooming::Bool

    dragstart::Vector{Float64}
    pointer::Vector{Float64}
    scale::Float64
end

function plot!(julia::JuliaView, f::Function, w::Integer, h::Integer)
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
    julia.plot = plt

    return plt
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

    view =
        is_mandel ?
        MandelView(0.0im, 4.0, plt, pressed, zooming, dragstart, pointer, scale) :
        JuliaView(0.0im, 4.0, c, plt, pressed, zooming, dragstart, pointer, scale)

    @guarded draw(canvas) do widget
        ctx = getgc(canvas)
        h = height(canvas)
        w = width(canvas)

        plt = plot!(view, f, w, h)
        surface = CairoImageSurface(plt)

        set_source_surface(ctx, surface, 0, 0)
        paint(ctx)
    end

    mouse_left = GtkGestureClick(canvas)
    signal_connect(mouse_left, "pressed") do controller, n_press, x, y
        view.pressed = true
        view.dragstart = [x, y]
        return
    end

    signal_connect(mouse_left, "released") do args...
        on_release(args..., view, f)
    end

    mouse_motion = GtkEventControllerMotion(canvas)
    signal_connect(mouse_motion, "motion") do args...
        on_motion(args..., view)
    end

    # mouse_scroll =
    #     GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL, canvas)

    # signal_connect(mouse_scroll, "scroll") do args...
    #     on_zoom(args..., view, f)
    # end

    show(win)

    if !isinteractive()
        signal_connect(win, :close_request) do widget
            GC.gc()
        end
        @async Gtk4.GLib.glib_main()
    end

    return view
end

end