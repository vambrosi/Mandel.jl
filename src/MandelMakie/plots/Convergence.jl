using GLMakie

function iterate(f, x, λ, ε, type)
    for i in 0:200
        if abs(x) < ε
            if type == :sa
                return i - log(2, -log(abs(x)))
            elseif type == :saε
                return i - log(2, log(sqrt(λ), abs(x)))
            elseif type == :a
                return i - log(abs(λ), abs(x))
            else
                return NaN
            end
        end
        x = f(x, λ)
    end
    return Inf
end

function plot_heights(f, xrange, λrange, εrange, correct)
    fig = Figure(size = (600, 650))

    ax = Axis(fig[1, 1])
    sg = SliderGrid(
        fig[2, 1],
        (label = "λ", range = λrange, format = "{:.7f}", startvalue = λrange[1]),
        (label = "ε", range = εrange, format = "{:.7f}", startvalue = εrange[1]),
    )

    sliderobservables = [s.value for s in sg.sliders]
    ys_a = lift(sliderobservables...) do λ, ε
        iterate.(f, xrange, λ, ε, :a)
    end

    ys_sa = lift(sliderobservables...) do λ, ε
        iterate.(f, xrange, λ, ε, correct ? :saε : :sa)
    end

    label = correct ? L"n - \log_2(\log_{\sqrt{\lambda}}(d))" : L"n - \log_2(-\log(d))"
    lines!(ax, xrange, ys_a, label = L"n - \log_\lambda (d))")
    lines!(ax, xrange, ys_sa, label = label)

    limits!(ax, xrange[1], xrange[end], -10, 10)
    axislegend(ax)

    return fig
end

let
    f(x, λ) = λ * x + x^2
    xrange = range(-1.0, 0.5, 10000)
    λrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    εrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    plot_heights(f, xrange, λrange, εrange, false)
end

let
    f(x, λ) = λ * x + x^2
    xrange = range(-1.0, 0.5, 10000)
    λrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    εrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    plot_heights(f, xrange, λrange, εrange, true)
end

let
    f(x, λ) = - λ * x + x^2
    xrange = range(-1.0, 0.5, 10000)
    λrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    εrange = 0.5 .* 10 .^ range(0.0, -7.0, 100)
    plot_heights(f, xrange, λrange, εrange, true)
end