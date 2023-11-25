# Plots for Central Limit Theorem proof
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using Random
using LaTeXStrings
using Plots

# plots a series
function plot_normal_mgf()
    t = range(start = -3, stop = 3, length = 1000)
    pl = plot(
        t,
        t -> exp(t^2 / 2),
        label = false,
        xlabel = L"t",
        ylabel = L"M_X(t)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    display(pl)
end


plot_normal_mgf()
savefig("~/Downloads/normal_mgf.svg")
