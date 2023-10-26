# Demonstrates the birthday paradox

# 2023 by Ralf Herbrich
# Hasso Plattner Institut

using Plots
using LaTeXStrings

M = 365
N = 80

function P(n, M)
    P = prod(range(M - n + 1, M) / M)
    return 1 - P
end

# plot the borthday paradox
p = plot(
    1:N,
    map(n -> P(n, M), 1:N),
    legend = false,
    linewidth = 3,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 16,
    yguidefontsize = 16,
)
scatter!([23], [P(23, M)])
ylabel!(L"P(\bar{A})")
xlabel!(L"n")
display(p)
