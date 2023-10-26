# Plots for Stirlings approximation
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Random
using LaTeXStrings
using Plots

# computes the Stirling approximation of log (n choose k)
function stirling(n, k)
    p = k / n
    H = -(p * log(p) + (1 - p) * log(1 - p))
    return (n * H)
end

# compute the true value of log (n choose k)
function log_n_choose_k(n, k)
    return log(factorial(n)) - log(factorial(k)) - log(factorial(n - k))
end

# plots the Stirling approximation of log (n choose k) for k = 1, ..., n-1
function plot_stirling_approximation(n)
    k = 1:(n-1)
    p = plot(
        k,
        [stirling(n, i) for i in k],
        label = "Stirling's Approximation",
        xlabel = L"k",
        ylabel = L"\log\binom{n}{k}",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    plot!(k, [log_n_choose_k(n, i) for i in k], linewidth = 3, label = "True Value")
    scatter!(k, [stirling(n, i) for i in k], color = :blue, markersize = 3)
    scatter!(
        k,
        [log_n_choose_k(n, i) for i in k],
        legend = false,
        color = :red,
        markersize = 3,
    )
    display(p)
end

function plot_proof_idea(n)
    k = 1:(n-1)
    p = plot(
        xlabel = L"k",
        ylabel = L"\log(k)",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        legend = false,
        xlims = (0, n),
    )
    for i in k
        x = [i, i, i + 1, i + 1]
        y = [0, log(i), log(i), 0]
        plot!(Shape(collect(zip(x, y))), label = false, color = :blue, alpha = 0.5)
    end
    k_hires = collect(1:0.1:n)
    y_hires = [log(i) for i in k_hires]
    k_hires_upper_bound = collect(1:0.1:n)
    y_hires_upper_bound = [log(i) for i in k_hires_upper_bound]
    k_hires_lower_bound = collect(2:0.1:n)
    y_hires_lower_bound = [log(i - 1) for i in k_hires_lower_bound]
    push!(k_hires, k_hires[end])
    push!(y_hires, 0)
    plot!(Shape(collect(zip(k_hires, y_hires))), alpha = 0.3)
    plot!(k_hires_upper_bound, y_hires_upper_bound, color = :red, linewidth = 3)
    plot!(k_hires_lower_bound, y_hires_lower_bound, color = :red, linewidth = 3)
    display(p)
end

plot_stirling_approximation(20)
plot_proof_idea(20)
