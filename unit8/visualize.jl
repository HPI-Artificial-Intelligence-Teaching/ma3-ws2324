# Plots for data visualization
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots
using CSV
using DataFrames

# plots a bar plot for a fixed dataset
function plot_starting_salary()
    salary = [57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 70]
    freq = [4, 1, 3, 5, 8, 10, 0, 5, 2, 3, 1]

    pl = bar(
        salary,
        freq,
        label = false,
        xlabel = L"\mathrm{Gehalt}",
        ylabel = L"\mathrm{Häufigkeit}",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# visualize the Old Faithful geyser dataset
function plot_histogram(filename = "/tmp/old_faithful.tsv"; no_bins = 8)
    # read the TSV file into a dataframe
    df = CSV.read(filename, delim = '\t', DataFrame)

    pl = histogram(
        df.eruptions,
        label = false,
        bins = no_bins,
        xlabel = L"\mathrm{Eruptionsdauer}",
        ylabel = L"\mathrm{Häufigkeit}",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# visualize the Old Faithful geyser dataset as a histogram video
function plot_histogram_video(filename = "/tmp/old_faithful.tsv")
    # read the TSV file into a dataframe
    df = CSV.read(filename, delim = '\t', DataFrame)
    (mn, mx) = extrema(df.eruptions)

    anim = @animate for bins = 8:300
        pl = histogram(
            df.eruptions,
            label = false,
            bins = range(start = mn, stop = mx, length = bins),
            xlabel = L"\mathrm{Eruptionsdauer}",
            ylabel = L"\mathrm{Häufigkeit}",
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
        )
        ylims!(0, 30)
        title!("Bins = $bins")
        display(pl)
    end

    return (anim)
end

plot_starting_salary()
savefig("~/Downloads/plot_starting_salary.svg")

plot_histogram(no_bins = 8)
savefig("~/Downloads/plot_histogram_8.svg")
plot_histogram(no_bins = 500)
savefig("~/Downloads/plot_histogram_500.svg")
anim = plot_histogram_video()
gif(anim, "~/Downloads/histogram.gif", fps = 30)
