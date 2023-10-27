using Zygote
using ChainRulesCore

@unpack σ, transf_matrices, u_fn = ts
@unpack global_dofs, solver = u_fn
@unpack penalty, problem, xmin = solver
@unpack Kes = solver.elementinfo
getdh(p::StiffnessTopOptProblem) = p.ch.dh

dh = getdh(problem)

σ = ts(PseudoDensities(x0))

n_dofs = length(global_dofs)

global_dofs
u_fn.x


f1 = Figure(resolution = (800,800))
ax1 = Axis(f1[1,1]) 

for i in keys(node_points)
    scatter!(ax1, node_points[i][1], node_points[i][2])
end

function f1(x::Vector{Float64})
    return x.^2 + 2*x .+ 1
end

function f2(x::Float64)
    return [x^4+2*x, x^2, 5*x]
end

using ChainRulesCore

function f3(x::Float64)
    return f1(f2(x))
end


function ChainRulesCore.rrule( f3::Function, y::Float64)
    val = f3(y)
    println("Hi from ChainRulesCore")
    function f1_pullback(Δ)
        return nothing, Tangent{typeof(f3)}(; y = Δ'*())
    end
    return val, f1_pullback
end

jacobian(f3,1.0)

