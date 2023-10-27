using ChainRulesCore
using ChainRulesTestUtils
using Revise
using UnPack
using TopOpt
using Ferrite
using Ferrite: getnbasefunctions, CellIterator

Zygote.gradient(ts)
test_rrule(ts, PseudoDensities(x0)[1,:]);


# ts::TrussStress{T}, x::PseudoDensities) where {T}
#or?
# ts::TrussStress, x::PseudoDensities)
@unpack σ, transf_matrices, u_fn = ts
@unpack global_dofs, solver = u_fn
@unpack penalty, problem, xmin = solver
@unpack Kes = solver.elementinfo
# K = ak(Kes)
dh = problem.ch.dh

# Forward-pass
σ = ts(PseudoDensities(x0))
u = u_fn(PseudoDensities(x0))
uu = u.u
# n_dofs = length(global_dofs)
getA(sp::TrussProblem) = [cs.A for cs in sp.truss_grid.crosssecs]
As = getA(problem)
# function truss_stress_pullback(Δ)
#       #gradient will be the same size as elements, or stress itself.
    Δσ = Vector{Float64}(undef, length(σ))
    Δu = Vector{Float64}(undef, 2*length(σ))
    # dudAs = jacobian(u, x)
    for e in 1:length(Δσ)
        celldofs!(global_dofs, dh, e)
        println(global_dofs)

        #σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
        # @show global_dofs
        # @show u[global_dofs]
        # @show transf_matrices[e]
        # @show (transf_matrices[e] * Kes[e]) / As[e]
        # Δσ[e] = norm((transf_matrices[e] * Kes[e]) / As[e] * duAs[global_dofs])
        println(size(Δσ))
        println("Transf matrix")
        @show size(transf_matrices[e])
        @show transf_matrices[e]
        println("Kes e")
        println(size(Kes[e]))
        @show Kes[e]

        @show Δu[global_dofs]
        @show Δu[global_dofs] = (transf_matrices[e] * Kes[e])[1,:]'/ As[e] * Δσ[e]

    end