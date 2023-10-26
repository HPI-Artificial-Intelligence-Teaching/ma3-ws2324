# Plots for histograms of uniform distributions
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

using Random
using LaTeXStrings
using Plots

# generate a vector of random 32-bit strings
function gen_bitarray(n)
    return [rand(Bool, 32) for i = 1:n]
end

# transform a bit array to a 32 bit integer
function bitarray_to_int(x)::UInt32
    return sum([x[i] * 2^(i - 1) for i in eachindex(x)])
end

# transform a bit array to a 32-bit floating point number
function bitarray_to_float(x)
    sign = (1 - 2 * x[32])
    exponent = bitarray_to_int(x[24:31]) - 127
    mantissa = 1 + bitarray_to_int(x[1:23]) / 2^23
    return sign * mantissa * 2.0^exponent
end


D = gen_bitarray(10^6)
p = histogram(
    map(bitarray_to_int, D),
    bins = 64,
    label = "integer number",
    xlabel = L"x",
    ylabel = L"P(x)",
    legend = :topleft,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 16,
    yguidefontsize = 16,
)
display(p)

p = histogram(
    map(bitarray_to_float, D),
    bins = 64,
    label = "floating point number",
    xlabel = L"x",
    ylabel = L"P(x)",
    legend = :topleft,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 16,
    yguidefontsize = 16,
)
display(p)

p = histogram(
    log10.(abs.(map(bitarray_to_float, D))),
    bins = 64,
    label = "floating point number",
    xlabel = L"\log_{10}(|x|)",
    ylabel = L"P\left(\log_{10}(|x|)\right)",
    legend = :topleft,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 16,
    yguidefontsize = 16,
)
display(p)
