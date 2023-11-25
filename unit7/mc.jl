# Plots for Monte Carlo integration
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# plots a series
function plot_quarter_circle()
    x = range(start = 0, stop = 1, length = 1000)
    pl = plot(
        x,
        x -> sqrt(1 - x^2),
        label = L"h(x)",
        xlabel = L"x",
        ylabel = L"y",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        aspect_ratio = :equal,
    )

    shp = [(0.0, 0.0)]
    append!(shp, collect(zip(x, map(t -> sqrt(1 - t^2), x))))
    plot!(Shape(shp), label = false, fillcolor = :blue, alpha = 0.2)
    xlims!(0, 1)
    ylims!(0, 1)
    display(pl)
end

# plots an evolving Monte-Carlo estiamte of a quarter area of a circle
function plot_monte_carlo_estimate(n; no_samples = 10, alpha = 0.2, color = :red)
    # plot the empty canvas
    p = plot(
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"n")
    ylabel!(L"\bar{h}(x)")

    # compute the indexing range
    ns = collect(1:n)

    # plot the samples
    for j = 1:no_samples
        xs = sqrt.(1 .- rand(Uniform(0, 1), n) .^ 2)
        h_bar = cumsum(xs) ./ ns
        plot!(ns, h_bar, linewidth = 1, alpha = alpha, color = color)
    end

    # plot the true value
    plot!(ns, n -> π / 4, linewidth = 2, color = :blue)


    display(p)
end


plot_quarter_circle()
savefig("~/Downloads/quarter_circle.svg")

plot_monte_carlo_estimate(1000)
savefig("~/Downloads/monte_carlo_estimate.svg")
