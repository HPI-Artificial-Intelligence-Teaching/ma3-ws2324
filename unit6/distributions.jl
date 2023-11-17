# Plots for continous distributions
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# plots a uniform distribution
function plot_uniform(a, b; xlims = nothing, max_y = 1.1)
    x = if (!isnothing(xlims))
        collect(xlims[1]:((xlims[2]-xlims[1])/1000):xlims[2])
    else
        collect((a-(b-a)/2):((b-a)/1000):(b+(b-a)/2))
    end
    d = Uniform(a, b)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(xlims))
        xlims!(xlims...)
    end
    ylims!(0, max_y)
    display(pl)
end

# plots an exponential distribution
function plot_exponential(λ; max_x = nothing)
    x = if (!isnothing(max_x))
        collect(0:max_x/1000:max_x)
    else
        collect(0:(2λ/1000):(2λ))
    end
    d = Exponential(λ)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(max_x))
        xlims!(0, max_x)
    end
    display(pl)
end

# plots the Euler-Gamma function
function plot_euler_gamma(; max_x = nothing, max_y = nothing)
    z = if (!isnothing(max_x))
        collect(0.15:max_x/1000:max_x)
    else
        collect(0.15:(6/1000):6)
    end
    pl = plot(
        z,
        z -> gamma(z),
        label = false,
        xlabel = L"z",
        ylabel = L"\Gamma(z)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(max_x))
        xlims!(0, max_x)
    end
    if (!isnothing(max_y))
        ylims!(0, max_y)
    end
    display(pl)
end

# plots a Gamma distribution
function plot_gamma(α, λ; max_x = nothing)
    x = if (!isnothing(max_x))
        collect(0:max_x/1000:max_x)
    else
        collect(0:(2λ/(α*1000)):(2λ))
    end
    d = Gamma(α, λ)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(max_x))
        xlims!(0, max_x)
    end
    display(pl)
end

# plot the Gamma distribution for different values of λ and fixed α
function plot_gamma_3D_fixed_alpha(α, λ_max)
    x = range(0, stop = 5, length = 100)
    λ = range(0.1, stop = λ_max, length = 100)
    x_coarse = range(0, stop = 5, length = 30)
    λ_coarse = range(0.1, stop = λ_max, length = 30)
    p = surface(
        x,
        λ,
        (x, λ) -> pdf(Gamma(α, λ), x),
        color = :reds,
        fillalpha = 0.1,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        ztickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        zguidefontsize = 16,
    )
    wireframe!(x_coarse, λ_coarse, (x, λ) -> pdf(Gamma(α, λ), x)), xlabel!(L"x")
    ylabel!(L"\lambda")
    zlabel!(L"p(x)")
    display(p)
end

# plot the Gamma distribution for different values of α and fixed λ
function plot_gamma_3D_fixed_lambda(λ, α_max)
    x = range(0, stop = 5, length = 100)
    α = range(1, stop = α_max, length = 100)
    x_coarse = range(0, stop = 5, length = 30)
    α_coarse = range(1, stop = α_max, length = 30)
    p = surface(
        x,
        α,
        (x, α) -> pdf(Gamma(α, λ), x),
        color = :reds,
        fillalpha = 0.1,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        ztickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        zguidefontsize = 16,
    )
    wireframe!(x_coarse, α_coarse, (x, α) -> pdf(Gamma(α, λ), x)), xlabel!(L"x")
    ylabel!(L"\alpha")
    zlabel!(L"p(x)")
    display(p)
end

# plots a Normal distribution
function plot_normal(μ, σ2; xlims = nothing, ylims = nothing)
    x = if (!isnothing(xlims))
        collect(xlims[1]:((xlims[2]-xlims[1])/1000):xlims[2])
    else
        collect((μ-3*sqrt(σ2)):((6*sqrt(σ2))/1000):(μ+3*sqrt(σ2)))
    end
    d = Normal(μ, sqrt(σ2))
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(xlims))
        xlims!(xlims...)
    end
    if (!isnothing(ylims))
        ylims!(ylims...)
    end
    display(pl)
end

# plots a Normal CDF distribution
function plot_normal_cdf(μ, σ2; xlims = nothing)
    x = if (!isnothing(xlims))
        collect(xlims[1]:((xlims[2]-xlims[1])/1000):xlims[2])
    else
        collect((μ-3*sqrt(σ2)):((6*sqrt(σ2))/1000):(μ+3*sqrt(σ2)))
    end
    d = Normal(μ, sqrt(σ2))
    pl = plot(
        x,
        x -> cdf(d, x),
        label = false,
        xlabel = L"t",
        ylabel = L"\Phi(t)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    if (!isnothing(xlims))
        xlims!(xlims...)
    end
    ylims!(0, 1)
    display(pl)
end

# plots a graphic to illustrate that Φ(-t)=1-Φ(t)
function plot_normal_symmetry(t)
    x = collect(-3:(6/1000):3)
    d = Normal()
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    low_tail_x = collect(range(start = -3, stop = t, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, t)
    push!(low_tail_x, -3)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )

    high_tail_x = collect(range(start = -t, stop = 3, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, 3)
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

# plots a 2D Normal distribution
function plot_2D_normal(
    μ,
    Σ;
    x1_lim = (-3, 3),
    x2_lim = (-3, 3),
    color = :greens,
    zlims = nothing,
)
    # plot the Normal distribution
    x1s = range(x1_lim[1], stop = x1_lim[2], length = 100)
    x2s = range(x2_lim[1], stop = x2_lim[2], length = 100)
    x1s_coarse = range(x1_lim[1], stop = x1_lim[2], length = 30)
    x2s_coarse = range(x2_lim[1], stop = x2_lim[2], length = 30)
    p = surface(
        x1s,
        x2s,
        (x1, x2) -> pdf(MvNormal(μ, Σ), [x1; x2])[1],
        color = color,
        fillalpha = 0.1,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        ztickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        zguidefontsize = 16,
    )
    wireframe!(x1s_coarse, x2s_coarse, (x1, x2) -> pdf(MvNormal(μ, Σ), [x1; x2])[1])
    xlabel!(L"x_1")
    ylabel!(L"x_2")
    zlabel!(L"p(x_1,x_2)")
    if (!isnothing(zlims))
        zlims!(zlims...)
    end
    display(p)
end

# plots the empiricial distribution of a sample os squared distances of standard normal deviations
function plot_empirical_chi_squared(n)
    d = Normal()
    d2 = Chisq(2)
    x = rand(d, n)
    y = rand(d, n)
    z = x .^ 2 + y .^ 2
    pl = histogram(
        z,
        label = false,
        normalize = :pdf,
        xlabel = L"X^2+Y^2",
        ylabel = L"\hat{p}(X^2+Y^2)",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    t = range(start = minimum(z), stop = maximum(z), length = 1000)
    plot!(t, t -> pdf(d2, t), linewidth = 3, color = :red)
    display(pl)
end

# plots a chi-squared distribution
function plot_chi_squared(n)
    d = Chisq(n)
    x = range(0, stop = 2 * n, length = 1000)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# plots a t-distribution
function plot_student_t(n; with_normal = false)
    d = TDist(n)
    d2 = Normal()
    x = range(-3, stop = 3, length = 1000)
    pl = plot(
        x,
        x -> pdf(d, x),
        label = false,
        xlabel = L"x",
        ylabel = L"p(x)",
        linewidth = 3,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    plot!(x, x -> pdf(d2, x), label = false, linewidth = 3, color = :red)

    display(pl)
end



plot_uniform(0, 1, xlims = (-0.5, 2.5), max_y = 1.1)
savefig("~/Downloads/uniform_0_1.svg")
plot_uniform(0, 2, xlims = (-0.5, 2.5), max_y = 1.1)
savefig("~/Downloads/uniform_0_2.svg")

plot_exponential(2, max_x = 5)
savefig("~/Downloads/exponential_2.svg")
plot_exponential(0.5, max_x = 5)
savefig("~/Downloads/exponential_0.5.svg")

plot_euler_gamma(max_x = 4.1, max_y = 8)
savefig("~/Downloads/euler_gamma.svg")

plot_gamma(0.5, 1, max_x = 5)
savefig("~/Downloads/gamma_0.5_1.svg")
plot_gamma(1, 1, max_x = 5)
savefig("~/Downloads/gamma_1_1.svg")
plot_gamma(2, 1, max_x = 5)
savefig("~/Downloads/gamma_2_1.svg")

plot_gamma(2, 0.5, max_x = 5)
savefig("~/Downloads/gamma_2_0.5.svg")
plot_gamma(2, 1, max_x = 5)
savefig("~/Downloads/gamma_2_1.svg")
plot_gamma(2, 1.5, max_x = 5)
savefig("~/Downloads/gamma_2_1.5.svg")

plot_gamma_3D_fixed_alpha(2, 1)
savefig("~/Downloads/gamma_fixed_alpha_2.png")
plot_gamma_3D_fixed_lambda(1, 3)
savefig("~/Downloads/gamma_fixed_lambda_1.png")

plot_normal(-1, 1, xlims = (-4, 4))
savefig("~/Downloads/normal_-1_1.svg")
plot_normal(0, 1, xlims = (-4, 4))
savefig("~/Downloads/normal_0_1.svg")
plot_normal(1, 1, xlims = (-4, 4))
savefig("~/Downloads/normal_1_1.svg")

plot_normal(0, 0.5, xlims = (-4, 4), ylims = (0, 0.7))
savefig("~/Downloads/normal_0_0.5.svg")
plot_normal(0, 1, xlims = (-4, 4), ylims = (0, 0.7))
savefig("~/Downloads/normal_0_1_v2.svg")
plot_normal(0, 2, xlims = (-4, 4), ylims = (0, 0.7))
savefig("~/Downloads/normal_0_2.svg")

plot_normal_cdf(0, 1, xlims = (-4, 4))
savefig("~/Downloads/normal_cdf_0_1.svg")
plot_normal_symmetry(-1.5)
savefig("~/Downloads/normal_symmetry.svg")
plot_2D_normal([0, 0], [1 0; 0 1])
savefig("~/Downloads/normal_2D_0_0_1_0_0_1_free.png")
plot_2D_normal([0, 0], [1 0; 0 1], zlims = (0, 0.38))
savefig("~/Downloads/normal_2D_0_0_1_0_0_1.png")
plot_2D_normal([0, 0], [2 0; 0 1], zlims = (0, 0.38))
savefig("~/Downloads/normal_2D_0_0_2_0_0_1.png")
plot_2D_normal([0, 0], [1 0.6; 0.6 1], zlims = (0, 0.38))
savefig("~/Downloads/normal_2D_0_0_1_0.6_0.6_1.png")
plot_2D_normal([0, 0], [1 0.9; 0.9 1], zlims = (0, 0.38))
savefig("~/Downloads/normal_2D_0_0_1_0.9_0.9_1.png")

plot_empirical_chi_squared(10000)
savefig("~/Downloads/empirical_chi_squared.svg")
plot_chi_squared(3)
savefig("~/Downloads/chi_squared_3.svg")
plot_chi_squared(10)
savefig("~/Downloads/chi_squared_10.svg")
plot_chi_squared(50)
savefig("~/Downloads/chi_squared_50.svg")

plot_student_t(1, with_normal = true)
savefig("~/Downloads/student_t_1.svg")
plot_student_t(10, with_normal = true)
savefig("~/Downloads/student_t_10.svg")
plot_student_t(100, with_normal = true)
savefig("~/Downloads/student_t_100.svg")
