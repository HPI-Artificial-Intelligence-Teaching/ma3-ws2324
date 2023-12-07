# Plots for parameter estimation
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# plots a histogram of (fictious) body heights as well as some estimators
function plot_body_heights(
    no_samples = 300;
    μ = 170,
    σ = 8,
    μσs = [(170, 8), (170, 9)],
    cols = [:red, :blue],
)
    # generate a random sample of body heights
    Random.seed!(1234)
    heights = rand(Normal(μ, σ), no_samples)

    pl = histogram(
        heights,
        bins = 20,
        normalize = true,
        label = false,
        xlabel = L"\mathrm{Größe}",
        ylabel = L"p(\mathrm{Größe})",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    for i in eachindex(μσs)
        (μ, σ) = μσs[i]
        col = cols[(i-1)%length(cols)+1]
        plot!(
            range(start = minimum(heights), stop = maximum(heights), length = 300),
            x -> pdf(Normal(μ, σ), x),
            label = false,
            linewidth = 3,
            color = col,
        )
    end

    μ_hat = mean(heights)
    σ_hat = std(heights)
    plot!(
        range(start = minimum(heights), stop = maximum(heights), length = 300),
        x -> pdf(Normal(μ_hat, σ_hat), x),
        label = false,
        linewidth = 3,
        color = :green,
    )

    println("Mean of the sample: ", μ_hat)
    println("Standard deviation of the sample: ", σ_hat)
    display(pl)
end

# plots the MSE of an estimator
function plot_estimator_mse(;
    θs = range(1, stop = 2, length = 100),
    dist = θ -> Normal(0, sqrt(θ)),
    no_samples = 10,
    estimators = [xs -> sum(xs .^ 2) / no_samples, xs -> sum(xs .^ 2) / (no_samples + 2)],
    n = 500000,
)
    # generate a random sample of body heights
    Random.seed!(1234)

    # compute the MSE for each θ
    function compute_mse(θ, est)
        mse = 0
        for i = 1:n
            xs = rand(dist(θ), no_samples)
            θ_hat = est(xs)
            mse += (θ_hat - θ)^2
        end
        mse /= n
    end

    # generate the plots on top of each other
    pl = plot(
        label = false,
        xlabel = L"\theta",
        ylabel = L"\mathrm{MSE}(\theta)",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    for estimator in estimators
        mse = map(θ -> compute_mse(θ, estimator), θs)
        plot!(θs, mse, label = false, linewidth = 3)
    end
    display(pl)
end

# plots the bias of an estimator
function plot_estimator_bias(;
    θs = range(1, stop = 2, length = 100),
    dist = θ -> Normal(0, sqrt(θ)),
    no_samples = 10,
    estimators = [xs -> sum((xs .- mean(xs)).^ 2) / no_samples, xs -> sum((xs .- mean(xs)).^ 2) / (no_samples - 1)],
    n = 1000000,
)
    # generate a random sample of body heights
    Random.seed!(1234)

    # compute the bias for each θ
    function compute_bias(θ, est)
        bias = 0
        for i = 1:n
            xs = rand(dist(θ), no_samples)
            θ_hat = est(xs)
            bias += θ_hat - θ
        end
        bias /= n
    end

    # generate the plots on top of each other
    pl = plot(
        label = false,
        xlabel = L"\theta",
        ylabel = L"E[T(X) - \theta)]",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    for estimator in estimators
        bias = map(θ -> compute_bias(θ, estimator), θs)
        plot!(θs, bias, label = false, linewidth = 3)
    end
    display(pl)
end


# plots the bias of an estimator
function plot_estimator_distribution(;
    dist = Normal(20, 4),
    no_samples = 10,
    estimators = [xs -> mean(xs), xs -> sum((xs .- mean(xs)).^ 2) / (no_samples - 1)],
    true_pdf = [x -> pdf(Normal(20, 4/sqrt(no_samples)),x), x -> (no_samples - 1) / 4^2 * pdf(Chisq(no_samples - 1), (no_samples - 1) / 4^2 * x)],
    file_names = ["~/Downloads/plot_estimator_distribution_mean.svg", "~/Downloads/plot_estimator_distribution_variance.svg"],
    n = 100000,
)
    # generate a random sample of body heights
    Random.seed!(1234)

    # generate the random samples
    X = Matrix{Float64}(undef, n, length(estimators))
    for i = 1:n
        x = rand(dist, no_samples)
        for j = 1:length(estimators)
            X[i,j] = estimators[j](x)
        end
    end

    # generate the plots 
    for j = 1:length(estimators)
        pl = histogram(
            X[:, j],
            label = false,
            normalize = true,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
        )
        plot!(
            range(start = minimum(X[:, j]), stop = maximum(X[:, j]), length = 300),
            x -> true_pdf[j](x),
            label = false,
            linewidth = 3,
            color = :red,
        )
        display(pl)
        savefig(file_names[j])
    end
end

# plots the likelihood for a sample
function plot_likelihood(x; dist = θ -> Normal(θ, 1), θs = range(-2, stop = 5, length = 100))
    pl = plot(
        θs,
        θ -> pdf(dist(θ), x),
        linewidth = 3,
        label = false,
        xlabel = L"\theta",
        ylabel = L"\mathcal{L}(\theta)",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end

# plots a random sample from a uniform distribution with a barplot
function plot_uniform_barplot(n; dist = Uniform(0, 3))
    Random.seed!(1234)
    xs = rand(dist, n)
    println("Random waiting time samples are: ", xs)
    println("MLE estimator: ", maximum(xs))
    println("MoM-1 estimator: ", 2*mean(xs))
    println("MoM-2 estimator: ", sqrt(3/n*sum(xs.^2)) )
    pl = bar(
        xs,
        fill(1/n, n),
        label = false,
        normalize = true,
        xlabel = L"x",
        ylabel = L"p(x)",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(pl)
end



plot_body_heights()
savefig("~/Downloads/plot_body_heights.svg")

plot_estimator_mse()
savefig("~/Downloads/plot_estimator_mse.svg")
plot_estimator_mse(
    θs = range(-2, stop = 2, length = 100),
    dist = θ -> Normal(θ, 1),
    no_samples = 10,
    estimators = [xs -> sum(xs) / 10, xs -> sum(xs) / 12],
    n = 500000,
)
savefig("~/Downloads/plot_estimator_mse2.svg")


plot_estimator_bias()
savefig("~/Downloads/plot_estimator_bias.svg")

plot_estimator_distribution()

plot_likelihood(1; dist = θ -> Normal(θ, 1), θs=range(-2, stop = 5, length = 1000))
savefig("~/Downloads/plot_likelihood_normal.svg")
plot_likelihood(1; dist = θ -> Uniform(0, θ), θs=range(0.1, stop = 5, length = 1000))
savefig("~/Downloads/plot_likelihood_uniform.svg")
plot_likelihood(1; dist = θ -> Uniform(θ, θ+1), θs=range(-2, stop = 5, length = 1000))
savefig("~/Downloads/plot_likelihood_uniform2.svg")

plot_uniform_barplot(5)
savefig("~/Downloads/plot_uniform_barplot.svg")