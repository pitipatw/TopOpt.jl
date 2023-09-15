using ChainRulesCore
using ChainRulesTestUtils
using UnPack
using TopOpt

function ChainRulesCore.rrule(
    ts::TrussStress, x::PseudoDensities)
    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    @unpack Kes = solver.elementinfo
    σ # u = ts(x)
    dh = getdh(problem)
    
    # Forward-pass
    truss_stress = ts(x)
    n_dofs = length(global_dofs)
    
    function truss_stress_pullback(Δ)
        #gradient will be the same size as elements, or stress itself.
        Δσ = Vector{Float64}(undef, length(σ))
        u = u_fn(x)
        for e in 1:length(σ)
            celldofs!(global_dofs, dh, e)
            #σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
            # Δσ[e] = Tangent{typeof(ts)}(;  -(transf_matrices[e] * Kes[e])[1] / As[e])

            _, dρe = get_ρ_dρ(x.x[e], penalty, xmin)
            celldofs!(global_dofs, dh, e)
            Keu = bcmatrix(Kes[e]) * u.u[global_dofs]
            dσdx_tmp[e] =  dot(Keu, solver.lhs[global_dofs])*u
        end 
        return nothing, Tangent{typeof(x)}(; dσdx_tmp)
    end
    return σ , truss_stress_pullback
end


test_rrule(ts, x0)