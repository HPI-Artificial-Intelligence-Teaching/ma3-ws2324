# Plots for Cauchy distribution
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Random
using LaTeXStrings
using Plots

# plots a Cauchy distribution
function plot_cauchy(n; min_x = -10, max_x = 10)
    x = range(min_x, max_x, length = n)
    p = plot(
        x,
        x -> 1 / π * 1 / (1 + x^2),
        label = false,
        xlabel = L"y",
        ylabel = L"f_Y(y)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(p)
end

plot_cauchy(1000)
