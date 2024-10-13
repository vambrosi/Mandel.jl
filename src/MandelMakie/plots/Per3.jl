using GLMakie

function height(c, R, N)
    z = 0

    for i in 1:N
        if abs(z) > R
            return i / 3 - log(2, log(abs(z)))
        end

        z = (z^2 - 1 - c + c^3) / (z^2 - c^2)
    end

    return NaN
end

function mod_iterate(c, R, N)
    z = 0

    for i in 1:N
        if abs(z) > R
            return mod(i, 3) * 0.5 / 3
        end

        z = (z^2 - 1 - c + c^3) / (z^2 - c^2)
    end

    return 0.5
end

function iteration_matrix(center, radius, pixels)
    pixel_length = 2 * radius / pixels
    corner = center - (radius + pixel_length / 2) * complex(1, 1)

    heights = Matrix{Float64}(undef, pixels, pixels)
    for j in 1:pixels
        for i in 1:pixels
            c = corner + pixel_length * complex(i, j)
            heights[i, j] = height(c, 1000.0, 1000) / 20.0

            isfinite(heights[i, j]) || (heights[i, j] = 0.5)
        end
    end

    return heights
end

function mod_matrix(center, radius, pixels)
    pixel_length = 2 * radius / pixels
    corner = center - (radius + pixel_length / 2) * complex(1, 1)

    heights = Matrix{Float64}(undef, pixels, pixels)
    for j in 1:pixels
        for i in 1:pixels
            c = corner + pixel_length * complex(i, j)
            heights[i, j] = mod_iterate(c, 1000.0, 1000)

            isfinite(heights[i, j]) || (heights[i, j] = 0.5)
        end
    end

    return heights
end

function Δ(f)
    Δf = zero(f)
    for y in 2:size(f, 2)-1, x in 2:size(f, 1)-1
        Δf[x, y] = f[x-1, y] + f[x+1, y] + f[x, y-1] + f[x, y+1] - 4f[x, y]
    end
    return Δf
end

function pool(A, a)
    B = zero(A)
    for y in 2:size(A, 2)-1
        for x in 2:size(A, 1)-1
            if maximum(abs, A[x-1:x+1, y-1:y+1] .- A[x, y]) > a
                B[x, y] = 0.5
            end
        end
    end
    return B
end

threshold(x, a) = x > a ? 0.5 : 0.0

function mask(A, t1, t2, t3, b, keep, all)
    B = zero(A)
    for y in 2:size(A, 2)-1
        for x in 2:size(A, 1)-1
            C = A[x-1:x+1, y-1:y+1]

            check =
                all ? (t1 ∈ C && t2 ∈ C && t3 ∈ C) :
                (t1 ∈ C && t2 ∈ C) || (t2 ∈ C && t3 ∈ C) || (t1 ∈ C && t3 ∈ C)
            
            if check
                B[x, y] = b
            else
                B[x, y] = keep ? A[x, y] : 0.0
            end
        end
    end
    return B
end

function plot_matrix(m)
    d = size(m, 1)

    fig = Figure(size = (d, d))
    ax = Axis(fig[1, 1], aspect = AxisAspect(1))

    heatmap!(ax, m, colormap = :twilight, colorrange = (0.0, 1.0))

    hidedecorations!(ax)
    hidespines!(ax)

    return fig
end

fig1 = plot_matrix(mod.(iteration_matrix(0.25, 3.0, 3000), 1.0));
save(joinpath("plots", "continuous.png"), fig1)

levels1 = threshold.(mod.(Δ(iteration_matrix(0.25, 3.0, 3000)), 1.0), 0.001);
fig2 = plot_matrix(levels1);
save(joinpath("plots", "laplacian.png"), fig2)

levels2 = pool(iteration_matrix(0.25, 3.0, 3000), 0.05);
fig3 = plot_matrix(levels2);
save(joinpath("plots", "pooling.png"), fig3)

fig4 = plot_matrix(mod_matrix(0.25, 3.0, 3000));
save(joinpath("plots", "preperiod_mod3.png"), fig4)

test_levels = [i * 0.5 / 3 for i in 0:2]
levels5 = mask(mod_matrix(0.25, 3.0, 3000), test_levels..., 0.75, true, false);
fig5 = plot_matrix(levels5);
save(joinpath("plots", "preperiod_mod3_border_at_least_2.png"), fig5)

levels6 = mask(mod_matrix(0.25, 3.0, 3000), test_levels..., 0.75, true, true);
fig6 = plot_matrix(levels6);
save(joinpath("plots", "preperiod_mod3_border_all_3.png"), fig6)

levels7 = mask(mod_matrix(0.25, 3.0, 3000), test_levels..., 0.75, false, false);
fig7 = plot_matrix(levels7);
save(joinpath("plots", "preperiod_mod3_only_border_at_least_2.png"), fig7)

levels8 = mask(mod_matrix(0.25, 3.0, 3000), test_levels..., 0.75, false, true);
fig8 = plot_matrix(levels8);
save(joinpath("plots", "preperiod_mod3_only_border_all_3.png"), fig8)