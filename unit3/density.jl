# Plots for uniform distribution 
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Random
using LaTeXStrings
using Plots

# plots a uniform distribution over the unit circle
function plot_uniform_distribution(n)
    ω = range(-π, 2.5π, length = n)
    p = plot(
        ω,
        ω -> (0 < ω ≤ 2π) ? 1 / 2π : 0,
        label = "Dichte über Winkel",
        xlabel = L"\omega",
        ylabel = L"p(\omega)",
        linewidth = 2,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    ω_full = [0; 0; 2π; 2π]
    p_full = [0; 1 / 2π; 1 / 2π; 0]
    plot!(
        Shape(collect(zip(ω_full, p_full))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    ω_partial = [π; π; 3 / 2 * π; 3 / 2 * π]
    p_partial = [0; 1 / 2π; 1 / 2π; 0]
    plot!(
        Shape(collect(zip(ω_partial, p_partial))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    display(p)
end

# plots a uniform distribution over the numbers 1 to n
function plot_discrete_uniform_cdf(n = 6)
    xs = range(-1, n + 1, length = 1000)
    p = plot(
        xs,
        x -> (0 < x ≤ n) ? round(x) / n : ((x < 0) ? 0 : 1),
        seriestype = :scatter,
        label = "Kumulative Verteilungsfunktion",
        xlabel = L"x",
        ylabel = L"F(x)",
        linewidth = 2,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(p)
end

# plots a uniform distribution over the unit circle (with CDF)
function plot_uniform_distribution_cdf(n)
    ω = range(-π, 2.5π, length = n)
    p = plot(
        ω,
        ω -> (0 < ω ≤ 2π) ? 1 / 2π : 0,
        label = "Dichte über Winkel",
        xlabel = L"\omega",
        ylabel = L"p(\omega), F(\omega)",
        linecolor = :orange,
        linewidth = 2,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    plot!(
        ω,
        ω -> (0 < ω ≤ 2π) ? ω / 2π : ((ω < 0) ? 0 : 1),
        linecolor = :blue,
        label = "Verteilungsfunktion über Winkel",
        xlabel = L"\omega",
        linewidth = 2,
    )
    display(p)
end

plot_uniform_distribution(1000)
plot_discrete_uniform_cdf(6)
plot_uniform_distribution_cdf(1000)
