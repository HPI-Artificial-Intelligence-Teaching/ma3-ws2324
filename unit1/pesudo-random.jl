# Plots for histograms of pseudo-random numbers
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

# struct that stores the state of a linear congruential generator
mutable struct LGR
    x::UInt32
end

# sets up a new linear congruential generator with a seed
function random(seed)
    return LGR(seed)
end

# implenents the C++ 11 pseudo-random number generator  
function rand(d::LGR)
    d.x = (16807 * d.x) % (2^31 - 1)
    return d.x
end

# generate 1,000,000 random numbers with the seed 42
using Plots

rnd = random(42)
D = [rand(rnd) for i = 1:10^6]
p = histogram(
    D ./ 2^31,
    label = :none,
    bins = 50,
    legend = :topleft,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 16,
    yguidefontsize = 16,
)
plot!([0, 1], [10^6 / 50, 10^6 / 50], seriestype = :line, linewidth = 3, label = :none)
display(p)

@time for i = 1:10^6
    rand(rnd)
end
