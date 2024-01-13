# Plots for regressions
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Distributions
using Random
using LaTeXStrings
using Plots
using Printf

# plots data together with linear fit
function plot_data(;
    x = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
    y = [45, 52, 54, 63, 62, 68, 75, 76, 92, 88],
    output_fit = true,
    output_summary = false,
    plot_abline = true,
    plot_residuals = false,
    residual_filename = nothing,
    filename = nothing,
)

    # produce the plot
    p = plot(
        x,
        y,
        seriestype = :scatter,
        xlabel = L"x",
        ylabel = L"y",
        label = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    # compute the fit
    S_XY = sum((x .- mean(x)) .* (y .- mean(y)))
    S_XX = sum((x .- mean(x)) .* (x .- mean(x)))
    S_YY = sum((y .- mean(y)) .* (y .- mean(y)))

    β_1 = S_XY / S_XX
    β_0 = mean(y) - β_1 * mean(x)

    # output the fit
    if output_fit
        @printf("Coefficients:\n%11s\t%8s\n%11.4f\t%8.4f\n", "(Intercept)", "x", β_0, β_1)
    end

    # output a summary
    if output_summary
        ŷ = β_0 .+ β_1 .* x
        SSR = sum((y .- ŷ) .* (y .- ŷ))
        σ̂² = SSR / (length(x) - 2)

        @printf("\nResiduals:\n")
        if (length(x) <= 7)
            for i in eachindex(x)
                @printf("%9d ", i)
            end
            @printf("\n")
            for i in eachindex(x)
                @printf("%9.6f ", y[i] - ŷ[i])
            end
        else
            @printf("%7s\t%7s\t%7s\t%6s\t%6s\n", "Min", "1Q", "Median", "3Q", "Max")
            @printf(
                "%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n",
                quantile(y .- ŷ, 0.00),
                quantile(y .- ŷ, 0.25),
                quantile(y .- ŷ, 0.50),
                quantile(y .- ŷ, 0.75),
                quantile(y .- ŷ, 1.00)
            )
        end

        @printf("\nCoefficients:\n             Estimate Std. Error t value Pr(>|t|)\n")
        @printf(
            "%-11s %9.5f %10.5f %7.3f %7.6f\n",
            "(Intercept)",
            β_0,
            sqrt((σ̂² * sum(x .^ 2)) / (length(x) * S_XX)),
            β_0 / sqrt((σ̂² * sum(x .^ 2)) / (length(x) * S_XX)),
            2 * (
                1 - cdf(
                    TDist(length(x) - 2),
                    abs(β_0 / sqrt((σ̂² * sum(x .^ 2)) / (length(x) * S_XX))),
                )
            )
        )
        @printf(
            "%-11s %9.5f %10.5f %7.3f %7.6f\n",
            "x",
            β_1,
            sqrt(σ̂² / S_XX),
            β_1 / sqrt(σ̂² / S_XX),
            2 * (1 - cdf(TDist(length(x) - 2), abs(β_1 / sqrt(σ̂² / S_XX))))
        )
        @printf("\n")
        @printf(
            "\nResidual standard error: %6.4f on %d degrees of freedom\n",
            sqrt(SSR / (length(x) - 2)),
            length(x) - 2
        )
        @printf("\n")
    end

    # plot the fit
    if plot_abline
        plot!(x, β_0 .+ β_1 .* x, label = false, linewidth = 2)
    end

    # plot the residuals
    if plot_residuals
        for i in eachindex(x)
            plot!(
                [x[i], x[i]],
                [y[i], β_0 + β_1 * x[i]],
                label = false,
                linewidth = 2,
                linecolor = :black,
            )
        end
    end

    display(p)
    if (!isnothing(filename))
        savefig(filename)
    end

    p2 = plot(
        x,
        y .- (β_0 .+ β_1 .* x),
        seriestype = :scatter,
        xlabel = L"x",
        ylabel = L"y - \hat{y}",
        label = false,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )
    display(p2)

    # if there is a filename, save the residual plot
    if (!isnothing(residual_filename))
        savefig(residual_filename)
    end
end

# plots the distribution of the MLE of the intercept, slope, and SSR
function plot_distribution_of_estimators(
    type;
    x = 1:15,
    β_0 = 2,
    β_1 = 3,
    σ = 2,
    n = 100000,
)
    S_XX = sum((x .- mean(x)) .* (x .- mean(x)))

    # simulate n estiamtions from zero-mean random noise
    estimate = Vector{Float64}(undef, n)
    for i = 1:n
        y = β_0 .+ β_1 .* x .+ rand(Normal(0, σ), length(x))
        S_XY = sum((x .- mean(x)) .* (y .- mean(y)))
        S_YY = sum((y .- mean(y)) .* (y .- mean(y)))
        if type == :intercept
            estimate[i] = mean(y) - (S_XY / S_XX) * mean(x)
        elseif type == :slope
            estimate[i] = S_XY / S_XX
        elseif type == :ssr
            estimate[i] = S_YY - (S_XY * S_XY) / S_XX
        end
    end

    # plot the histogram 
    xlabel_string =
        (type == :intercept) ? L"$\hat{\beta}_0$" :
        ((type == :slope) ? L"$\hat{\beta}_1$" : L"$\mathrm{SSR}$")
    ylabel_string =
        (type == :intercept) ? L"$\hat{p}\left(\hat{\beta}_0\right)$" :
        (
            (type == :slope) ? L"$\hat{p}\left(\hat{\beta}_1\right)$" :
            L"$\hat{p}\left(\mathrm{SSR}\right)$"
        )

    p = histogram(
        estimate,
        label = false,
        normalize = :pdf,
        xlabel = xlabel_string,
        ylabel = ylabel_string,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
    )

    # add the true distributions
    f = if type == :intercept
        t -> pdf(Normal(β_0, sqrt((σ^2 * sum(x .^ 2)) / (length(x) * S_XX))), t)
    elseif type == :slope
        t -> pdf(Normal(β_1, σ / √S_XX), t)
    elseif type == :ssr
        t -> 1 / σ^2 * pdf(Chisq(length(x) - 2), t / σ^2)
    end

    plot!(
        range(start = minimum(estimate), stop = maximum(estimate), length = 300),
        x -> f(x),
        linewidth = 3,
        color = :red,
        legend = false,
    )

    display(p)
end

plot_data(
    x = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
    y = [45, 52, 54, 63, 62, 68, 75, 76, 92, 88],
    filename = "~/Downloads/regress_artificial_data.svg",
)
plot_data(
    x = [1, 2, 3, 4, 5, 6, 7],
    y = [3, 2, 5, 6, 4, 8, 9],
    plot_residuals = true,
    filename = "~/Downloads/regress_fiber.svg",
)
plot_data(
    x = [46, 53, 29, 61, 36, 39, 47, 49, 52, 38, 55, 32, 57, 54, 44],
    y = [12, 15, 7, 17, 10, 11, 11, 12, 14, 9, 16, 8, 18, 14, 12],
    filename = "~/Downloads/regress_moisture.svg",
)
plot_data(
    x = [45, 50, 55, 60, 65, 70, 75],
    y = [24.2, 25, 23.3, 22, 21.5, 20.6, 19.8],
    output_summary = true,
    filename = "~/Downloads/regress_consumption.svg",
)
plot_data(
    x = [60, 62, 64, 65, 66, 67, 68, 70, 72, 74],
    y = [63.6, 65.2, 66, 65.5, 66.9, 67.1, 67.4, 68.3, 70.1, 70],
    output_summary = true,
    residual_filename = "~/Downloads/regress_father_son_residuals.svg",
    filename = "~/Downloads/regress_father_son.svg",
)

Random.seed!(42)
x = rand(Uniform(0, 20), 30)
y = x .^ 2 .+ rand(Normal(0, 5), length(x))
plot_data(
    x = x,
    y = y,
    residual_filename = "~/Downloads/regress_quadratic_residuals.svg",
    filename = "~/Downloads/regress_quadratic.svg",
)
plot_data(
    x = x .^ 2,
    y = y,
    residual_filename = "~/Downloads/regress_input_transform_residuals.svg",
    filename = "~/Downloads/regress_input_transform_features.svg",
)
y = 10 .* exp.(x ./ 4) .+ rand(Normal(0, 3), length(x))
plot_data(
    x = x,
    y = y,
    residual_filename = "~/Downloads/regress_exponential_residuals.svg",
    filename = "~/Downloads/regress_exponential.svg",
)
plot_data(
    x = x,
    y = log.(y),
    residual_filename = "~/Downloads/regress_exponential_output_transform_residuals.svg",
    filename = "~/Downloads/regress_exponential_output_transform.svg",
)


plot_distribution_of_estimators(:intercept)
savefig("~/Downloads/regress_intercept.svg")
plot_distribution_of_estimators(:slope)
savefig("~/Downloads/regress_slope.svg")
plot_distribution_of_estimators(:ssr)
savefig("~/Downloads/regress_ssr.svg")
