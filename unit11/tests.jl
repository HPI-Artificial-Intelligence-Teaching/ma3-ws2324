# Plots for hypothesis tests
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using DataFrames
using CSV
using SpecialFunctions
using Random
using LaTeXStrings
using Plots

# read the ImDB files
function read_imdb_data(
    title_url = "https://datasets.imdbws.com/title.basics.tsv.gz",
    ratings_url = "https://datasets.imdbws.com/title.ratings.tsv.gz",
)
    # download the files
    println("Downloading files...")
    basics_file = download(title_url)
    ratings_file = download(ratings_url)

    # read the files
    println("Reading files...")
    basics = CSV.read(basics_file, delim = '\t', quoted = false, DataFrame)
    ratings = CSV.read(ratings_file, delim = '\t', quoted = false, DataFrame)

    # join the files
    println("Joining files...")
    titles = innerjoin(basics, ratings, on = :tconst)

    # return the data frame
    titles
end

# plots the histogram for the ratings 
function plot_ratings_histogram(titles; normal_fit = true, bins = :sqrt)
    p = histogram(
        titles.averageRating,
        xlabel = "Rating",
        ylabel = L"\hat{p}(\mathrm{Rating})",
        bins = bins,
        label = false,
        normalize = :pdf,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    if (normal_fit)
        μ = mean(titles.averageRating)
        σ = std(titles.averageRating)
        plot!(
            range(
                minimum(titles.averageRating),
                maximum(titles.averageRating),
                length = 1000,
            ),
            x -> pdf(Normal(μ, σ), x),
            linewidth = 3,
            color = :red,
            label = false,
        )
    end

    display(p)
end

# helper function to get the indistinguishable movies
function get_indistinguishable(table; α = 0.05)
    μ = mean(table.averageRating)
    σ = √var(table.averageRating)
    n = length(table.averageRating)
    z = quantile(Normal(0, 1), α / 2)
    μs = μ .+ [z * σ / √n, -z * σ / √n]
    println("μ = ", μ)
    println("σ = ", σ)
    println("n = ", n)
    println("μs = ", μs)
    return (filter(row -> row.averageRating > μs[1] && row.averageRating < μs[2], table))
end

# plots the errors of a hypothesis tests
function plot_errors_gauss_test(;
    μ_range = range(5, 6, length = 300),
    σ = √3,
    n = 266,
    α = 0.05,
    μ0 = 5.6,
)
    # compute the accept region
    z_lower = μ0 + σ / √n * quantile(Normal(0, 1), α / 2)
    z_upper = μ0 + σ / √n * quantile(Normal(0, 1), 1 - α / 2)

    # plot the power
    p = plot(
        μ_range,
        μ -> cdf(Normal(μ, σ / √n), z_lower) + (1 - cdf(Normal(μ, σ / √n), z_upper)),
        xlabel = L"\mathrm{Wahres\ } \mu",
        ylabel = L"P(\mathrm{Fehler})",
        linewidth = 3,
        label = "Typ-I Fehler",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    plot!(
        μ_range,
        μ -> cdf(Normal(μ, σ / √n), z_upper) - cdf(Normal(μ, σ / √n), z_lower),
        linewidth = 3,
        color = :red,
        label = "Typ-II Fehler",
    )

    plot!(
        [minimum(μ_range), maximum(μ_range)],
        [α, α],
        linewidth = 1,
        linestyle = :dash,
        color = :black,
        label = "Irrtumsniveau",
    )
    ylims!(0, 1)

    display(p)
end

# plots an inverted two-sided interval range
function plot_inverted_two_sided_interval_normal_distributions(; μ = 0, σ = 1, α = 0.05)
    # temporary variables
    d = Normal(μ, σ)
    z = quantile(Normal(0, 1), 0.5 + α / 2)

    # generate the feasiable x-region
    xs = range(start = μ - 3σ, stop = μ + 3σ, length = 1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color = :black,
        linewidth = 3,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    xlabel!(L"x")
    ylabel!(L"\mathcal{N}(x;0,1)")

    # plots the tails
    low_tail_x = collect(range(start = μ - 3σ, stop = μ - z * σ, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, μ - z * σ)
    push!(low_tail_x, μ - 3σ)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )

    high_tail_x = collect(range(start = μ + z * σ, stop = μ + 3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ + 3σ)
    push!(high_tail_x, μ + z * σ)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = μ - z * σ, stop = μ + z * σ, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, μ + z * σ)
    push!(mid_tail_x, μ - z * σ)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )

    display(p)
end

# plots the errors of a weak hypothesis tests
function plot_errors_weak_gauss_test(;
    μ_range = range(5, 6, length = 300),
    σ = √3,
    n = 266,
    α = 0.05,
    μ0 = 5.6,
)
    # compute the reject region
    z_lower = μ0 + σ / √n * quantile(Normal(0, 1), 0.5 - α / 2)
    z_upper = μ0 + σ / √n * quantile(Normal(0, 1), 0.5 + α / 2)

    # plot the power
    p = plot(
        μ_range,
        μ -> cdf(Normal(μ, σ / √n), z_upper) - cdf(Normal(μ, σ / √n), z_lower),
        xlabel = L"\mathrm{Wahres\ } \mu",
        ylabel = L"P(\mathrm{Fehler})",
        linewidth = 3,
        label = "Typ-I Fehler",
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    plot!(
        μ_range,
        μ -> cdf(Normal(μ, σ / √n), z_lower) + (1 - cdf(Normal(μ, σ / √n), z_upper)),
        linewidth = 3,
        color = :red,
        label = "Typ-II Fehler",
    )

    plot!(
        [minimum(μ_range), maximum(μ_range)],
        [α, α],
        linewidth = 1,
        linestyle = :dash,
        color = :black,
        label = "Irrtumsniveau",
    )
    ylims!(0, 1)

    display(p)
end

# plots the distribution of the Wilcoxon rank-sum statistic
function plot_wilcoxon_rank_sum_distribution(; n1 = 5, n2 = 5, approx = true)
    # recursively computes the distribution of the Wilcoxon rank-sum statistic
    function wilcoxon_rank_sum_distribution(k; n1 = 5, n2 = 5)
        if (k < 0 || n1 < 0 || n2 < 0)
            return 0
        end

        if (n1 == 1 && n2 == 0)
            return (k <= 0) ? 0 : 1
        end

        if (n1 == 0 && n2 == 1)
            return (k < 0) ? 0 : 1
        end

        P1 = wilcoxon_rank_sum_distribution(k - (n1 + n2), n1 = n1 - 1, n2 = n2)
        P2 = wilcoxon_rank_sum_distribution(k, n1 = n1, n2 = n2 - 1)

        return n1 / (n1 + n2) * P1 + n2 / (n1 + n2) * P2
    end

    k = 0:((n1+1)*n1/2+n1*n2)
    P = [wilcoxon_rank_sum_distribution(ki, n1 = n1, n2 = n2) for ki in k]
    p = diff([0; P])

    # plot the probability mass function
    pl = plot(
        k,
        p,
        xlabel = L"W",
        ylabel = L"p(W)",
        seriestype = :bar,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    if (approx)
        # plot the approximation
        plot!(
            k,
            k -> pdf(Normal(n1 * (n1 + n2 + 1) / 2, √(n1 * n2 * (n1 + n2 + 1) / 12)), k),
            linewidth = 3,
            color = :red,
            label = false,
        )
    end

    display(pl)
end

# computes the W statistic for two vectors
function compute_W(xs, ys)
    U = 0
    n = length(xs)
    for x in xs
        for y in ys
            U += (x > y) ? 1 : ((x == y) ? 0.5 : 0)
        end
    end
    return U + n * (n + 1) / 2
end

# plots a two-sided interval range for a chi-squared distribution
function plot_two_sided_interval_χ2_distributions(; n = 10, α = 0.05)
    # temporary variables
    d = Chisq(n - 1)
    z_lower = quantile(d, α / 2)
    z_upper = quantile(d, 1 - α / 2)
    μ = n - 1
    σ = sqrt(2 * (n - 1))

    # generate the feasiable x-region
    xs = range(start = 0, stop = μ + 3σ, length = 1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color = :black,
        linewidth = 3,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 12,
        yguidefontsize = 12,
    )
    xlabel!(L"\lambda_{\mathrm{LR}}")
    ylabel!(L"\mathcal{\chi^2}(\lambda_{\mathrm{LR}})")

    # plots the tails
    low_tail_x = collect(range(start = 0, stop = z_lower, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, z_lower)
    push!(low_tail_x, 0)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [z_lower / 2, z_lower / 2],
        [4 * pdf(d, z_lower / 2), 0.5 * pdf(d, z_lower / 2)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        z_lower / 2,
        5 * pdf(d, z_lower / 2),
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    high_tail_x = collect(range(start = z_upper, stop = μ + 3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ + 3σ)
    push!(high_tail_x, z_upper)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [(z_upper + μ + 3σ) / 2, (z_upper + μ + 3σ) / 2],
        [4 * pdf(d, (z_upper + μ + 3σ) / 2), 0.5 * pdf(d, (z_upper + μ + 3σ) / 2)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        (z_upper + μ + 3σ) / 2,
        5 * pdf(d, (z_upper + μ + 3σ) / 2),
        L"\frac{\alpha}{2}",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = z_lower, stop = z_upper, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, z_upper)
    push!(mid_tail_x, z_lower)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ + σ, μ],
        [0.9 * pdf(d, μ), 0.8 * pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        μ + σ,
        pdf(d, μ),
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots a lower one-sided interval range for a chi-squared distribution
function plot_lower_one_sided_interval_χ2_distributions(; n = 10, α = 0.05)
    # temporary variables
    d = Chisq(n - 1)
    z_lower = quantile(d, α)
    μ = n - 1
    σ = sqrt(2 * (n - 1))

    # generate the feasiable x-region
    xs = range(start = 0, stop = μ + 3σ, length = 1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color = :black,
        linewidth = 3,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 12,
        yguidefontsize = 12,
    )
    xlabel!(L"\lambda_{\mathrm{LR}}")
    ylabel!(L"\mathcal{\chi^2}(\lambda_{\mathrm{LR}})")

    # plots the tails
    low_tail_x = collect(range(start = 0, stop = z_lower, length = 300))
    low_tail_y = pdf.(d, low_tail_x)
    push!(low_tail_x, z_lower)
    push!(low_tail_x, 0)
    push!(low_tail_y, 0)
    push!(low_tail_y, 0)
    plot!(
        Shape(collect(zip(low_tail_x, low_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [z_lower / 2, z_lower / 2],
        [4 * pdf(d, z_lower / 2), 0.5 * pdf(d, z_lower / 2)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        z_lower / 2,
        5 * pdf(d, z_lower / 2),
        L"\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = z_lower, stop = μ + 3σ, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, μ + 3σ)
    push!(mid_tail_x, z_lower)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ + σ, μ],
        [0.9 * pdf(d, μ), 0.8 * pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        μ + σ,
        pdf(d, μ),
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# plots an upper one-sided interval range for a chi-squared distribution
function plot_upper_one_sided_interval_χ2_distributions(; n = 10, α = 0.05)
    # temporary variables
    d = Chisq(n - 1)
    z_upper = quantile(d, 1 - α)
    μ = n - 1
    σ = sqrt(2 * (n - 1))

    # generate the feasiable x-region
    xs = range(start = 0, stop = μ + 3σ, length = 1000)

    # generate the normal distribution plot
    p = plot(
        xs,
        x -> pdf(d, x),
        color = :black,
        linewidth = 3,
        legend = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 12,
        yguidefontsize = 12,
    )
    xlabel!(L"\lambda_{\mathrm{LR}}")
    ylabel!(L"\mathcal{\chi^2}(\lambda_{\mathrm{LR}})")

    high_tail_x = collect(range(start = z_upper, stop = μ + 3σ, length = 300))
    high_tail_y = pdf.(d, high_tail_x)
    push!(high_tail_x, μ + 3σ)
    push!(high_tail_x, z_upper)
    push!(high_tail_y, 0)
    push!(high_tail_y, 0)
    plot!(
        Shape(collect(zip(high_tail_x, high_tail_y))),
        label = false,
        fillcolor = :blue,
        alpha = 0.5,
    )
    plot!(
        [(z_upper + μ + 3σ) / 2, (z_upper + μ + 3σ) / 2],
        [4 * pdf(d, (z_upper + μ + 3σ) / 2), 0.5 * pdf(d, (z_upper + μ + 3σ) / 2)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        (z_upper + μ + 3σ) / 2,
        5 * pdf(d, (z_upper + μ + 3σ) / 2),
        L"\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    # plots the high-density region
    mid_tail_x = collect(range(start = 0, stop = z_upper, length = 300))
    mid_tail_y = pdf.(d, mid_tail_x)
    push!(mid_tail_x, z_upper)
    push!(mid_tail_x, 0)
    push!(mid_tail_y, 0)
    push!(mid_tail_y, 0)
    plot!(
        Shape(collect(zip(mid_tail_x, mid_tail_y))),
        label = false,
        fillcolor = :red,
        alpha = 0.5,
    )
    plot!(
        [μ + σ, μ],
        [0.9 * pdf(d, μ), 0.8 * pdf(d, μ)],
        color = :black,
        linewidth = 1,
        arrow = arrow(:closed),
    )
    annotate!(
        μ + σ,
        pdf(d, μ),
        L"1-\alpha",
        halign = :left,
        valign = :center,
        color = :black,
        rotation = 0,
        fontsize = 12,
    )

    display(p)
end

# performs a likelihood ratio test for two Normal samples
function lr_test(x, y; σ = √3)
    μ = mean([x; y])
    logL = sum(log.(pdf.(Normal(μ, σ), [x; y])))
    println("μ = ", μ, ",  logL = ", logL)

    μ_x = mean(x)
    logL_x = sum(log.(pdf.(Normal(μ_x, σ), x)))
    println("μ_x = ", μ_x, ",  logL = ", logL_x)
    μ_y = mean(y)
    logL_y = sum(log.(pdf.(Normal(μ_y, σ), y)))
    println("μ_y = ", μ_y, ",  logL = ", logL_y)

    λ = -2 * (logL - (logL_x + logL_y))
    return λ
end

# simulate the distribution of p-values
function plot_pvalue_distribution(; type = :both, n = 1000000)
    # generate the data
    samples = rand(Normal(0, 1), n)

    if (type == :both)
        # plot the histogram
        p = histogram(
            2 .* cdf.(Normal(0, 1), -abs.(samples)),
            xlabel = L"$\alpha(X)$",
            ylabel = L"$\hat{p}(\alpha(X))$",
            label = false,
            normalize = :pdf,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
        )
    elseif (type == :left)
        # plot the histogram
        p = histogram(
            cdf.(Normal(0, 1), samples),
            xlabel = L"$\alpha(X)$",
            ylabel = L"$\hat{p}(\alpha(X))$",
            label = false,
            normalize = :pdf,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
        )
    elseif (type == :right)
        # plot the histogram
        p = histogram(
            1 .- cdf.(Normal(0, 1), samples),
            xlabel = L"$\alpha(X)$",
            ylabel = L"$\hat{p}(\alpha(X))$",
            label = false,
            normalize = :pdf,
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
        )
    end
    plot!(
        range(0, 1, length = 1000),
        x -> pdf(Uniform(0, 1), x),
        linewidth = 3,
        color = :red,
        label = false,
    )

    display(p)
end



# uncomment this line to download the data
titles = read_imdb_data()

scifi2022_movies = filter(
    row ->
        contains(row.genres, "Sci-Fi") &&
            tryparse(Int64, row.startYear) == 2022 &&
            row.titleType == "movie",
    titles,
)
plot_ratings_histogram(scifi2022_movies, bins = :sqrt)
savefig("~/Downloads/scifi2022_movies_histogram.svg")
println(filter(row -> row.primaryTitle == "Jurassic World Dominion", scifi2022_movies))
get_indistinguishable(scifi2022_movies)

scifi2021_movies = filter(
    row ->
        contains(row.genres, "Sci-Fi") &&
            tryparse(Int64, row.startYear) == 2021 &&
            row.titleType == "movie",
    titles,
)
plot_ratings_histogram(scifi2021_movies, bins = :sqrt)
savefig("~/Downloads/scifi2021_movies_histogram.svg")

scifi_movies =
    filter(row -> contains(row.genres, "Sci-Fi") && row.titleType == "movie", titles)
plot_ratings_histogram(scifi_movies, bins = :scott)
savefig("~/Downloads/scifi_movies_histogram.svg")

romance2022_movies = filter(
    row ->
        contains(row.genres, "Romance") &&
            tryparse(Int64, row.startYear) == 2022 &&
            row.titleType == "movie",
    titles,
)
plot_ratings_histogram(romance2022_movies, bins = :auto)
savefig("~/Downloads/romance2022_movies_histogram.svg")
println(filter(row -> row.primaryTitle == "Marry Me", romance2022_movies))

x = mean(romance2022_movies.averageRating)
σ = std(romance2022_movies.averageRating, corrected = true)
n = length(romance2022_movies.averageRating)
stat =
    (
        x -
        filter(row -> row.primaryTitle == "Marry Me", romance2022_movies)[
            1,
            :,
        ].averageRating
    ) / (σ / √n)

println("x = ", x)
println("σ = ", σ)
println("n = ", n)
println("T = ", stat)
println(
    "two-sided: [",
    quantile(TDist(n - 1), 0.025),
    ", ",
    quantile(TDist(n - 1), 0.975),
    "]",
)
println("one-sided: [", quantile(TDist(n - 1), 0.05), ", +infty)")
println("one-sided: (-infty, ", quantile(TDist(n - 1), 0.95), "]")

plot_errors_gauss_test(μ_range = range(4.5, 6.5, length = 300))
savefig("~/Downloads/errors_gauss_test.svg")

plot_inverted_two_sided_interval_normal_distributions()
savefig("~/Downloads/inverted_two_sided_interval_normal_distributions.svg")

plot_errors_weak_gauss_test()
savefig("~/Downloads/errors_weak_gauss_test.svg")

plot_wilcoxon_rank_sum_distribution(n1 = 2, n2 = 1)
savefig("~/Downloads/wilcoxon_rank_sum_distribution_2_1.svg")
plot_wilcoxon_rank_sum_distribution(n1 = 3, n2 = 2)
savefig("~/Downloads/wilcoxon_rank_sum_distribution_3_2.svg")
plot_wilcoxon_rank_sum_distribution(n1 = 4, n2 = 4)
savefig("~/Downloads/wilcoxon_rank_sum_distribution_4_4.svg")
plot_wilcoxon_rank_sum_distribution(n1 = 8, n2 = 8)
savefig("~/Downloads/wilcoxon_rank_sum_distribution_8_8.svg")

n = length(scifi2021_movies.averageRating)
m = length(scifi2022_movies.averageRating)
W = compute_W(scifi2021_movies.averageRating, scifi2022_movies.averageRating)
println("n = ", n)
println("m = ", m)
println("W = ", W)
println("W_nm = ", √3 * (2 * W - n * (n + m + 1)) / √(n * m * (n + m + 1)))

plot_two_sided_interval_χ2_distributions(α = 0.05, n = 10)
savefig("~/Downloads/two_sided_interval_chi_squared_distributions.svg")
plot_lower_one_sided_interval_χ2_distributions(α = 0.05, n = 10)
savefig("~/Downloads/lower_one_sided_interval_chi_squared_distributions.svg")
plot_upper_one_sided_interval_χ2_distributions(α = 0.05, n = 10)
savefig("~/Downloads/upper_one_sided_interval_chi_squared_distributions.svg")

λ = lr_test(scifi2021_movies.averageRating, scifi2022_movies.averageRating)
println("λ = ", λ)
quantile(Chisq(1), 0.95)

plot_pvalue_distribution(type = :both)
savefig("~/Downloads/pvalue_distribution_both.svg")
plot_pvalue_distribution(type = :left)
savefig("~/Downloads/pvalue_distribution_left.svg")
plot_pvalue_distribution(type = :right)
savefig("~/Downloads/pvalue_distribution_right.svg")
