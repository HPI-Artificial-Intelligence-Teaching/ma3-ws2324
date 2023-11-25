# Plots for convergence
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# plots a series
function plot_series(f; n = 20, ϵ = nothing)
    x = collect(1:n)
    pl = plot(
        x,
        x -> f(x),
        label = false,
        xlabel = L"n",
        ylabel = L"x_n",
        seriestype = :scatter,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(ϵ))
        plot!(pl, x, x -> ϵ, label = false, color = :red, linewidth = 3)
    end
    display(pl)
end

# plots a graphic to illustrate that Φ(-t)=1-Φ(t)
function plot_normal_symmetry(t; μ = 0, σ = 1)
    x = collect(-4:(8/1000):4)
    d = Normal(μ, σ)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"z",
        ylabel = L"p(z)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    low_tail_x = collect(range(start = -4, stop = t, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, t)
    push!(low_tail_x, -4)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )

    high_tail_x = collect(range(start = -t, stop = 4, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, 4)
    push!(high_tail_x, -t)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )

    display(pl)
end

# plots a graphic to almost sure convergence
function plot_almost_sure_convergence(N)
    x = collect(range(start = -1.2, stop = 1.2, length = 1000))
    pl = plot(
        x,
        x -> (x < 0) ? 0 : 1,
        label = false,
        xlabel = L"x",
        ylabel = L"F(x)",
        linewidth = 3,
        color = :red,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    for n = 1:N
        d = Uniform(0, 1 / n)
        plot!(x, x -> cdf(d, x), label = false, linewidth = 1, color = :blue)
    end
    display(pl)
end

# plots a video to almost sure convergence
function plot_almost_sure_convergence_video(N)
    x = collect(range(start = -1.2, stop = 1.2, length = 1000))
    anim = @animate for n = 1:N
        pl = plot(
            x,
            x -> (x < 0) ? 0 : 1,
            label = L"X",
            xlabel = L"x",
            ylabel = L"F(x)",
            linewidth = 3,
            color = :red,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
            legendfontsize = 16,
        )

        d = Uniform(0, 1 / n)
        plot!(x, x -> cdf(d, x), label = L"X_{%$n}", linewidth = 3, color = :blue)
        display(pl)
    end

    return (anim)
end

# plots sample from a random process of a converging mean
function plot_mean_process(
    n;
    no_samples = 10,
    highlight_idx = 1,
    alpha = 0.1,
    color = :red,
    dist = Uniform(1.9, 2.1),
)
    # plot the empty canvas
    p = plot(
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"n")
    ylabel!(L"\bar{X}_n")

    # compute the indexing range
    ns = collect(1:n)

    for j = 1:no_samples
        width, α = if (j == highlight_idx)
            3, 1
        else
            1, alpha
        end
        xs = rand(dist, n)
        ys = cumsum(xs) ./ ns
        plot!(ns, ys, linewidth = width, alpha = α, color = color)
    end
    display(p)
end

# plots sample from a random process of a converging mean
function plot_mean_process_with_distribution_function(
    n;
    no_samples = 10,
    highlight_idx = 1,
    alpha = 0.1,
    color = :red,
    dist_function = (t, x) -> Uniform(1.9, 2.1),
)
    # plot the empty canvas
    p = plot(
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"n")
    ylabel!(L"\bar{X}_n")

    # compute the indexing range
    ns = collect(1:n)


    xs = Vector{Float64}(undef, n)
    for j = 1:no_samples
        width, α = if (j == highlight_idx)
            3, 1
        else
            1, alpha
        end
        old_x = 0.0
        for i = 1:n
            xs[i] = rand(dist_function(i, old_x), 1)[1]
            old_x = xs[i]
        end
        ys = cumsum(xs) ./ ns
        plot!(ns, ys, linewidth = width, alpha = α, color = color)
    end
    display(p)
end


# plots the empirical distribution of the mean
function plot_mean_distribution(n; no_samples = 5000, dist = Uniform(1.9, 2.1))
    # compute the samples of the means
    x_bar = Vector{Float64}(undef, no_samples)
    for j = 1:no_samples
        xs = rand(dist, n)
        x_bar[j] = mean(xs)
    end

    # plots the histogram
    p = histogram(
        x_bar,
        label = :none,
        legend = :none,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"Häufigkeit")
    xlabel!(L"\bar{X}_{%$n}")

    display(p)
end

# plots the empirical distribution of the mean with a distribution
function plot_mean_distribution_with_distribution_function(
    n;
    no_samples = 5000,
    dist_function = (t, x) -> Uniform(1.9, 2.1),
)
    # compute the samples of the means
    x_bar = Vector{Float64}(undef, no_samples)
    xs = Vector{Float64}(undef, n)
    for j = 1:no_samples
        old_x = 0.0
        for i = 1:n
            xs[i] = rand(dist_function(i, old_x), 1)[1]
            old_x = xs[i]
        end
        x_bar[j] = mean(xs)
    end

    # plots the histogram
    p = histogram(
        x_bar,
        label = :none,
        legend = :none,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"Häufigkeit")
    xlabel!(L"\bar{X}_{%$n}")

    display(p)
end




plot_series(n -> 1 / n, n = 20, ϵ = 0.2)
savefig("~/Downloads/series_1_n.svg")

plot_normal_symmetry(-1, μ = 0, σ = √2)
savefig("~/Downloads/normal_symmetry_-1_0_2.svg")

plot_almost_sure_convergence(10)
savefig("~/Downloads/almost_sure_convergence_10.svg")

anim = plot_almost_sure_convergence_video(100)
gif(anim, "~/Downloads/almost_sure_convergence.gif", fps = 30)

plot_mean_process(100, alpha = 0.2, no_samples = 10)
savefig("~/Downloads/mean_process_100.svg")
plot_mean_process(1000, alpha = 0.1, no_samples = 50)
savefig("~/Downloads/mean_process_1000.svg")
plot_mean_process(10000, alpha = 0.05, no_samples = 100)
savefig("~/Downloads/mean_process_10000.png")

plot_mean_process_with_distribution_function(
    10000,
    alpha = 0.2,
    no_samples = 10,
    dist_function = (j, x) -> Uniform(-j, j),
)
savefig("~/Downloads/mean_process_with_distribution_function_10000.svg")
plot_mean_process_with_distribution_function(
    10000,
    alpha = 0.1,
    no_samples = 10,
    dist_function = (j, x) -> Normal(x, 1),
)
savefig("~/Downloads/mean_process_with_distribution_function_10000_normal.svg")
plot_mean_process_with_distribution_function(
    10000,
    alpha = 0.1,
    no_samples = 1000,
    dist_function = (j, x) -> Normal(x, 1),
)
savefig("~/Downloads/mean_process_with_distribution_function_10000_normal.png")

plot_mean_distribution(2)
savefig("~/Downloads/mean_distribution_2.svg")
plot_mean_distribution(5)
savefig("~/Downloads/mean_distribution_5.svg")
plot_mean_distribution(10)
savefig("~/Downloads/mean_distribution_10.svg")
plot_mean_distribution(100)
savefig("~/Downloads/mean_distribution_100.svg")

plot_mean_distribution_with_distribution_function(
    1000,
    dist_function = (j, x) -> Normal(x, 10 / j),
)
savefig("~/Downloads/mean_distribution_with_distribution_function_1000_normal.svg")
plot_mean_distribution_with_distribution_function(
    10000,
    dist_function = (j, x) -> Normal(x, 10 / j),
)
savefig("~/Downloads/mean_distribution_with_distribution_function_10000_normal.svg")
plot_mean_distribution_with_distribution_function(
    1000,
    dist_function = (j, x) -> Uniform(-j, j),
)
savefig("~/Downloads/mean_distribution_with_distribution_function_1000_uniform.svg")
plot_mean_distribution_with_distribution_function(
    10000,
    dist_function = (j, x) -> Uniform(-j, j),
)
savefig("~/Downloads/mean_distribution_with_distribution_function_10000_uniform.svg")
