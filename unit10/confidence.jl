# Plots for confidence intervals
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# plots a series of normal distributions with different means
function plot_normal_distributions_across_means(;
    μs=range(start=-3, stop=3, length=100),
    σ=1,
    μss=[-2, -1, 0, 1, 2],
)
    # generate the feasiable x-region
    xs = range(start=minimum(μs) - 3σ, stop=maximum(μs) + 3σ, length=100)

    # generate a 3D plot of the sampling distribution
    μs_coarse = range(start=minimum(μs), stop=maximum(μs), length=30)
    xs_coarse = range(start=minimum(xs), stop=maximum(xs), length=30)
    p = surface(
        xs,
        μs,
        (x, μ) -> pdf(Normal(μ, σ), x),
        color=:greens,
        fillalpha=0.1,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        ztickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
        zguidefontsize=16,
        camera=(15, 30),
    )
    wireframe!(
        xs_coarse,
        μs_coarse, 
        (x, μ) -> pdf(Normal(μ, σ), x),
    )
    xlabel!(L"x")
    ylabel!(L"\mu")
    zlabel!(L"\mathcal{N}(x;\mu,\sigma^2)")

    # plot particular distributions
    for μ in μss
        μs = μ * ones(100)
        ps = map((x, μ) -> pdf(Normal(μ, σ), x), xs, μs)
        plot!(xs, μs, ps, color = :red, linewidth = 3)
    end

    display(p)
end

# plots a simple two-sided interval range
function plot_two_sided_interval_normal_distributions(;
    μ=0,
    σ=1,
    α=0.05,
)
    # temporary variables
    d = Normal(μ, σ)
    z = quantile(Normal(0, 1), 1 - α / 2)

    # generate the feasiable x-region
    xs = range(start=μ - 3σ, stop=μ + 3σ, length=1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color=:black,
        linewidth=3,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )
    xlabel!(L"x")
    ylabel!(L"\mathcal{N}(x;\mu,\sigma^2)")

    # plots the tails
    low_tail_x = collect(range(start = μ-3σ, stop = μ-z*σ, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, μ-z*σ)
    push!(low_tail_x, μ-3σ)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [μ-(3+z)/2*σ, μ-(3+z)/2*σ],
        [0.5*pdf(d, μ-(3+z)/2*σ)+0.1, 0.5*pdf(d, μ-(3+z)/2*σ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ-(3+z)/2*σ,
        0.5*pdf(d, μ-(3+z)/2*σ)+0.13,
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    high_tail_x = collect(range(start = μ+z*σ, stop = μ+3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ+3σ)
    push!(high_tail_x, μ+z*σ)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [μ+(3+z)/2*σ, μ+(3+z)/2*σ],
        [0.5*pdf(d, μ+(3+z)/2*σ)+0.1, 0.5*pdf(d, μ+(3+z)/2*σ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+(3+z)/2*σ,
        0.5*pdf(d, μ+(3+z)/2*σ)+0.13,
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = μ-z*σ, stop = μ+z*σ, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, μ+z*σ)
    push!(mid_tail_x, μ-z*σ)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ+σ, μ],
        [0.75*pdf(d, μ)+0.1, 0.75*pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+σ,
        0.75*pdf(d, μ)+0.11,
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots a simple two-sided confidence interval
function plot_two_sided_normal_distributions_confidence(;
    μs=range(start=-3, stop=3, length=100),
    σ=1,
    α=0.05,
    x = 0,
)
    # temporary variables
    z = quantile(Normal(0, 1), α / 2)
    μ_min = minimum(μs)
    μ_max = maximum(μs)

    # plot the tail regions
    low_tail_x = [μ_min-3σ, μ_min+z*σ, μ_max+z*σ, μ_min-3σ]
    low_tail_y = [μ_min, μ_min, μ_max, μ_max]
    p = plot(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )

    high_tail_x = [μ_min-z*σ, μ_max+3σ, μ_max+3σ, μ_max-z*σ]
    high_tail_y = [μ_min, μ_min, μ_max, μ_max]
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        fillcolor = :blue,
        alpha = 0.5,
    )

    # plot the high-density region
    mid_tail_x = [μ_min+z*σ, μ_min-z*σ, μ_max-z*σ, μ_max+z*σ]
    mid_tail_y = [μ_min, μ_min, μ_max, μ_max]
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        fillcolor = :red,
        alpha = 0.5,
    )

    # plot the confidence interval
    plot!(
        [x, x],
        [x + z*σ, x - z*σ],
        color = :black,
        linewidth = 3,
        arrow=arrow(:closed,:both),
    )
    println("z = ", z)
    println("x + z*σ = ", x + z*σ)
    println("x - z*σ = ", x - z*σ)

    # plot the labels
    xlabel!(L"x")
    ylabel!(L"\mu")

    display(p)
end

# plots a lower one-sided interval range
function plot_lower_one_sided_interval_normal_distributions(;
    μ=0,
    σ=1,
    α=0.05,
)
    # temporary variables
    d = Normal(μ, σ)
    z = quantile(Normal(0,1), α)
    
    # generate the feasiable x-region
    xs = range(start=μ - 3σ, stop=μ + 3σ, length=1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color=:black,
        linewidth=3,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )
    xlabel!(L"x")
    ylabel!(L"\mathcal{N}(x;\mu,\sigma^2)")

    # plots the tail
    low_tail_x = collect(range(start = μ-3σ, stop = μ+z*σ, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, μ+z*σ)
    push!(low_tail_x, μ-3σ)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [μ-(3-z)/2*σ, μ-(3-z)/2*σ],
        [0.5*pdf(d, μ-(3-z)/2*σ)+0.1, 0.5*pdf(d, μ-(3-z)/2*σ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ-(3-z)/2*σ,
        0.5*pdf(d, μ-(3-z)/2*σ)+0.13,
        L"\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = μ+z*σ, stop = μ+3*σ, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, μ+3*σ)
    push!(mid_tail_x, μ+z*σ)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ+σ, μ],
        [0.75*pdf(d, μ)+0.1, 0.75*pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+σ,
        0.75*pdf(d, μ)+0.11,
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots a lower one-sided confidence interval
function plot_lower_one_sided_normal_distributions_confidence(;
    μs=range(start=-3, stop=3, length=100),
    σ=1,
    α=0.05,
    x = 0,
)
    # temporary variables
    z = quantile(Normal(0, 1), α)
    μ_min = minimum(μs)
    μ_max = maximum(μs)

    # plot the tail regions
    low_tail_x = [μ_min-3σ, μ_min+z*σ, μ_max+z*σ, μ_min-3σ]
    low_tail_y = [μ_min, μ_min, μ_max, μ_max]
    p = plot(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )

    # plot the high-density region
    mid_tail_x = [μ_min+z*σ, μ_max+3σ, μ_max+3σ, μ_max+z*σ]
    mid_tail_y = [μ_min, μ_min, μ_max, μ_max]
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        fillcolor = :red,
        alpha = 0.5,
    )

    # plot the confidence interval
    plot!(
        [x, x],
        [μ_min, x-z*σ],
        color = :black,
        linewidth = 3,
        arrow=arrow(:closed),
    )

    # plot the labels
    xlabel!(L"x")
    ylabel!(L"\mu")

    display(p)
end

# plots an upper one-sided interval range
function plot_upper_one_sided_interval_normal_distributions(;
    μ=0,
    σ=1,
    α=0.05,
)
    # temporary variables
    d = Normal(μ, σ)
    z = quantile(Normal(0,1), 1 - α)
    
    # generate the feasiable x-region
    xs = range(start=μ - 3σ, stop=μ + 3σ, length=1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color=:black,
        linewidth=3,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )
    xlabel!(L"x")
    ylabel!(L"\mathcal{N}(x;\mu,\sigma^2)")

    # plots the tail
    high_tail_x = collect(range(start = μ+z*σ, stop = μ+3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ+3σ)
    push!(high_tail_x, μ+z*σ)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [μ+(3+z)/2*σ, μ+(3+z)/2*σ],
        [0.5*pdf(d, μ+(3+z)/2*σ)+0.1, 0.5*pdf(d, μ+(3+z)/2*σ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+(3+z)/2*σ,
        0.5*pdf(d, μ+(3+z)/2*σ)+0.13,
        L"\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = μ-3*σ, stop = μ+z*σ, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, μ+z*σ)
    push!(mid_tail_x, μ-3*σ)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ+σ, μ],
        [0.75*pdf(d, μ)+0.1, 0.75*pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+σ,
        0.75*pdf(d, μ)+0.11,
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots an upper one-sided confidence interval
function plot_upper_one_sided_normal_distributions_confidence(;
    μs=range(start=-3, stop=3, length=100),
    σ=1,
    α=0.05,
    x = 0,
)
    # temporary variables
    z = quantile(Normal(0, 1), 1 - α)
    μ_min = minimum(μs)
    μ_max = maximum(μs)

    # plot the tail regions
    low_tail_x = [μ_min+z*σ, μ_max+3*σ, μ_max+3*σ, μ_max+z*σ]
    low_tail_y = [μ_min, μ_min, μ_max, μ_max]
    p = plot(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )

    # plot the high-density region
    mid_tail_x = [μ_min-3σ, μ_min+z*σ, μ_max+z*σ, μ_min-3*σ]
    mid_tail_y = [μ_min, μ_min, μ_max, μ_max]
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        fillcolor = :red,
        alpha = 0.5,
    )

    # plot the confidence interval
    plot!(
        [x, x],
        [μ_max, x - z*σ],
        color = :black,
        linewidth = 3,
        arrow=arrow(:closed),
    )

    # plot the labels
    xlabel!(L"x")
    ylabel!(L"\mu")

    display(p)
end

# plots a two-sided interval range for a chi-squared distribution
function plot_two_sided_interval_χ2_distributions(;
    n=10,
    α=0.05,
)
    # temporary variables
    d = Chisq(n - 1)
    z_lower = quantile(d, α / 2)
    z_upper = quantile(d, 1 - α / 2)
    μ = n - 1
    σ = sqrt(2 * (n - 1))

    # generate the feasiable x-region
    xs = range(start=0, stop=μ+3σ, length=1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color=:black,
        linewidth=3,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=12,
        yguidefontsize=12,
    )
    xlabel!(L"(n-1) \cdot \hat{V}[X_1,\ldots,X_n] / \sigma^2")
    ylabel!(L"\mathcal{\chi^2}((n-1) \cdot \hat{V}[X_1,\ldots,X_n] / \sigma^2;n-1)")

    # plots the tails
    low_tail_x = collect(range(start = 0, stop = z_lower, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, z_lower)
    push!(low_tail_x, 0)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [z_lower/2, z_lower/2],
        [4*pdf(d, z_lower/2), 0.5*pdf(d, z_lower/2)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        z_lower/2,
        5 * pdf(d, z_lower/2),
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    high_tail_x = collect(range(start = z_upper, stop = μ+3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ+3σ)
    push!(high_tail_x, z_upper)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [(z_upper + μ + 3σ)/2, (z_upper + μ + 3σ)/2],
        [4*pdf(d, (z_upper + μ + 3σ)/2), 0.5*pdf(d, (z_upper + μ + 3σ)/2)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        (z_upper + μ + 3σ)/2,
        5 * pdf(d, (z_upper + μ + 3σ)/2),
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = z_lower, stop = z_upper, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, z_upper)
    push!(mid_tail_x, z_lower)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ+σ, μ],
        [0.9*pdf(d, μ), 0.8*pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow=arrow(:closed),
    )
    annotate!(
        μ+σ,
        pdf(d, μ),
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots a two-sided confidence interval for a chi-squared distribution
function plot_two_sided_χ2_distributions_confidence(;
    σs=range(start=0, stop=1, length=100),
    n=50,
    α=0.05,
    x = 0.5,
)
    # temporary variables
    z_lower = quantile(Chisq(n - 1), α / 2)
    z_upper = quantile(Chisq(n - 1), 1 - α / 2)
    σ_min = minimum(σs)
    σ_max = maximum(σs)
    x_max = σ_max*z_upper/(n-1)

    # plot the tail regions
    low_tail_x = [0, σ_min*z_lower/(n-1), σ_max*z_lower/(n-1), 0]
    low_tail_y = [σ_min, σ_min, σ_max, σ_max]
    p = plot(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )

    high_tail_x = [σ_min*z_upper/(n-1), x_max, x_max, σ_max*z_upper/(n-1)]
    high_tail_y = [σ_min, σ_min, σ_max, σ_max]
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        fillcolor = :blue,
        alpha = 0.5,
    )

    # plot the high-density region
    mid_tail_x = [σ_min*z_lower/(n-1), σ_min*z_upper/(n-1), σ_max*z_upper/(n-1), σ_max*z_lower/(n-1)]
    mid_tail_y = [σ_min, σ_min, σ_max, σ_max]
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        fillcolor = :red,
        alpha = 0.5,
    )

    # plot the confidence interval
    plot!(
        [x, x],
        [x * (n-1) / z_upper, x * (n-1) /  z_lower],
        color = :black,
        linewidth = 3,
        arrow=arrow(:closed,:both),
    )

    # plot the labels
    xlabel!(L"\hat{V}[X_1,\ldots,X_n]")
    ylabel!(L"\sigma^2")

    display(p)
end

# plot a binomial distribution
function plot_binomial_distribution(;
    n=30,
    p=0.1,
)
    d = Binomial(n, p)

    # generate the feasiable x-region
    xs = collect(0:n)

    # generate the normal distribution plot
    pl = plot(
        xs,
        x -> pdf(d, x),
        seriestype = :bar,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )

    # plot normal approximation
    d = Normal(n*p, sqrt(n*p*(1-p)))
    plot!(
        range(start = 0, stop = n, length = 300),
        x -> pdf(d, x),
        color=:red,
        linewidth=3,
    )
    xlabel!(L"y")
    ylabel!(L"\mathrm{Bin}(y;n,p)")

    display(pl)
end

# plot an approximate binomial confidence interval
function plot_approximate_binomial_confidence(;
    n = 50,
    α = 0.05,
    approx = false,
)
    # temporary variables
    z_upper = quantile(Normal(0, 1), 1 - α / 2)
    z_lower = quantile(Normal(0, 1), α / 2)

    # compute the right function for the confidence interval
    f = if (approx)
        (k, p) -> (k - n * p) / sqrt(k * (1 - k/n))
    else
        (k, p) -> (k - n * p) / sqrt(n * p * (1 - p))
    end

    # plot the Ω×Θ region
    ps = range(start = 0, stop = 1, length = 10*(n+1))
    ks = collect(0:n)
    pl = heatmap(
        ks,
        ps,
        (k, p) -> f(k,p) < z_upper && f(k,p) > z_lower ? 1 : 0,
        color=palette([:blue, :red], 2),
        fillalpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )
    xlabel!(L"k")
    ylabel!(L"p")

    display(pl)
end

# plot an approximate binomial confidence interval using Tschebyscheff's inequality
function plot_tschebycheff_approximate_binomial_confidence(;
    n = 50,
    α = 0.05,
    approx = false,
)
    # temporary variables
    z_upper = sqrt(1 / α)
    z_lower = -sqrt(1 / α)

    # compute the right function for the confidence interval
    f = (k, p) -> (k - n * p) / sqrt(n * p * (1 - p))

    # plot the Ω×Θ region
    ps = range(start = 0, stop = 1, length = 10*(n+1))
    ks = collect(0:n)
    pl = heatmap(
        ks,
        ps,
        (k, p) -> f(k,p) < z_upper && f(k,p) > z_lower ? 1 : 0,
        color=palette([:blue, :red], 2),
        fillalpha = 0.5,
        legend=false,
        xtickfontsize=14,
        ytickfontsize=14,
        xguidefontsize=16,
        yguidefontsize=16,
    )
    xlabel!(L"k")
    ylabel!(L"p")

    display(pl)
end


# checks that the confidence procedure is calibrated
function check_confidence(; n = 100000, α = 0.05)
    z = quantile(Normal(0, 1), 1 - α / 2)
    mistakes = 0
    for i = 1:n
        μ = rand(Uniform(-3,3))
        x = rand(Normal(μ, 1))
        if μ < x - z || μ > x + z
            mistakes += 1
        end    
    end
    println("Mistake rate: ", mistakes / n)
end

plot_normal_distributions_across_means()
savefig("~/Downloads/normal_distributions_across_means.png")
plot_two_sided_interval_normal_distributions(α=0.05, σ=1)
savefig("~/Downloads/two_sided_interval_normal_distributions.svg")
plot_two_sided_normal_distributions_confidence(α=0.05, σ=1, x=0.5)
savefig("~/Downloads/two_sided_normal_distributions_confidence.svg")
plot_lower_one_sided_interval_normal_distributions(α=0.05, σ=1)
savefig("~/Downloads/lower_one_sided_interval_normal_distributions.svg")
plot_lower_one_sided_normal_distributions_confidence(α=0.05, σ=1, x=0.5)
savefig("~/Downloads/lower_one_sided_normal_distributions_confidence.svg")
plot_upper_one_sided_interval_normal_distributions(α=0.05, σ=1)
savefig("~/Downloads/upper_one_sided_interval_normal_distributions.svg")
plot_upper_one_sided_normal_distributions_confidence(α=0.05, σ=1, x=0.5)
savefig("~/Downloads/upper_one_sided_normal_distributions_confidence.svg")
plot_two_sided_interval_χ2_distributions(α=0.05, n=10)
savefig("~/Downloads/two_sided_interval_chi_squared_distributions.svg")
plot_two_sided_χ2_distributions_confidence(α=0.05, n=50, x=0.64)
savefig("~/Downloads/two_sided_chi_squared_distributions_confidence.svg")
plot_binomial_distribution(p = 0.1)
savefig("~/Downloads/binomial_distribution_p=0.1.svg")
plot_binomial_distribution(p = 0.2)
savefig("~/Downloads/binomial_distribution_p=0.2.svg")

plot_approximate_binomial_confidence(n = 5, α = 0.05, approx = false)
savefig("~/Downloads/approximate_binomial_confidence_n=5_α=0.05_no_approximation.png")
plot_approximate_binomial_confidence(n = 5, α = 0.05, approx = true)
savefig("~/Downloads/approximate_binomial_confidence_n=5_α=0.05_approximation.png")

plot_approximate_binomial_confidence(n = 50, α = 0.05, approx = false)
savefig("~/Downloads/approximate_binomial_confidence_n=50_α=0.05_no_approximation.png")
plot_approximate_binomial_confidence(n = 50, α = 0.05, approx = true)
savefig("~/Downloads/approximate_binomial_confidence_n=50_α=0.05_approximation.png")

plot_approximate_binomial_confidence(n = 200, α = 0.05, approx = false)
savefig("~/Downloads/approximate_binomial_confidence_n=200_α=0.05_no_approximation.png")
plot_approximate_binomial_confidence(n = 200, α = 0.05, approx = true)
savefig("~/Downloads/approximate_binomial_confidence_n=200_α=0.05_approximation.png")

plot_tschebycheff_approximate_binomial_confidence(n = 50, α = 0.05)
savefig("~/Downloads/tschebycheff_approximate_binomial_confidence_n=50_α=0.05.png")

check_confidence()