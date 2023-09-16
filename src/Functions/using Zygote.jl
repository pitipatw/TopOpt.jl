using Zygote
using ChainRulesCore


function Y(x::Float64)
    return x^2
end

function Z(y::Float64)
    out = Vector{Float64}(undef,3)
    for i in eachindex(out)
        out[i] = y^i
    end
    return out
end

x = 2.0

Y(x)

p = Z(Y(x))

gradient(p)