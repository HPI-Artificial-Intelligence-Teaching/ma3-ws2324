# Plots for discrete distributions
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using Random
using LaTeXStrings
using Plots

# plots a binomial distribution
function plot_binomial(n, p)
    k = collect(0:n)
    d = Binomial(n, p)
    pl = plot(
        k,
        k -> pdf(d, k),
        label = false,
        xlabel = L"k",
        ylabel = L"p(k)",
        seriestype = :bar,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# plots a hyper-geometric distribution
function plot_hypergeometric(M, N, n)
    k = collect(0:min(M, n))
    d = Hypergeometric(M, N, n)
    pl = plot(
        k,
        k -> pdf(d, k),
        label = false,
        xlabel = L"k",
        ylabel = L"p(k)",
        seriestype = :bar,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# plots a Poisson distribution
function plot_poisson(λ; max_k = 10)
    k = collect(0:max_k)
    d = Poisson(λ)
    pl = plot(
        k,
        k -> pdf(d, k),
        label = false,
        xlabel = L"k",
        ylabel = L"p(k)",
        seriestype = :bar,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end


plot_binomial(10, 0.5)
savefig("~/Downloads/binomial_10_0.5.svg")
plot_binomial(10, 0.3)
savefig("~/Downloads/binomial_10_0.3.svg")
plot_binomial(10, 0.6)
savefig("~/Downloads/binomial_10_0.6.svg")

plot_hypergeometric(10, 10, 10)
savefig("~/Downloads/hypergeometric_10_10_10.svg")
plot_hypergeometric(100, 100, 10)
savefig("~/Downloads/hypergeometric_100_100_10.svg")
plot_hypergeometric(1000, 1000, 10)
savefig("~/Downloads/hypergeometric_1000_1000_10.svg")
plot_hypergeometric(6, 43, 6)
savefig("~/Downloads/hypergeometric_lotto.svg")

plot_poisson(4)
savefig("~/Downloads/poisson_4.svg")
plot_poisson(0.8)
savefig("~/Downloads/poisson_0.8.svg")
