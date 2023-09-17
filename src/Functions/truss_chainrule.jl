# TODO complete
# """
# rrule for autodiff.
    
# du/dxe = -K^-1 * dK/dxe * u
# d(u)/d(x_e) = - K^-1 * d(K)/d(x_e) * u
#             = - K^-1 * (Σ_ei d(ρ_ei)/d(x_e) * K_ei) * u
#             = - K^-1 * [d(ρ_e)/d(x_e) * K_e * u]
# d(u)/d(x_e)' * Δ = -d(ρ_e)/d(x_e) * u' * K_e * (K^-1 * Δ)
# """
function ChainRulesCore.rrule(
    ::typeof(ts), x::PseudoDensities)
    #or?
    # ts::TrussStress, x::PseudoDensities)
    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    @unpack Kes = solver.elementinfo
    K = ak(Kes)
    dh = getdh(problem)
    
    # Forward-pass
    σ = ts(x)
    u = u_fn(x)
    n_dofs = length(global_dofs)
    function truss_stress_pullback(Δ)
#         #gradient will be the same size as elements, or stress itself.
        Δσ = Vector{Float64}(undef, length(σ))
        for e in 1:length(Δσ)
            celldofs!(global_dofs, dh, e)
            #σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
            Δσ[e] = (transf_matrices[e] * Kes[e])[1] / As[e] * u[global_dofs]
        end
#             _, dρe = get_ρ_dρ(x.x[e], penalty, xmin)
#             celldofs!(global_dofs, dh, e)
#             Keu = bcmatrix(Kes[e]) * u.u[global_dofs]
#             dσdx_tmp[e] = -dρe * dot(Keu, solver.lhs[global_dofs])
#         end 

    # σ::AbstractVector{T} # stress vector, axial stress per cell
    # u_fn::Displacement
    # transf_matrices::AbstractVector{<:AbstractMatrix{T}}
    # fevals::Int
    # maxfevals::Int
    return Tangent{typeof(ts)}(;
        σ = Δσ,
        u_fn = u ,
        transf_matrices = NoTangent(),
        fevals = NoTangent(),
        maxfevals = NoTangent()), NoTangent()
end
    return σ , truss_stress_pullback
end

    # return truss_stress, Δ -> begin # v
    #     solver.rhs .= Δ
    #     solver(reuse_chol = true, assemble_f = false)
    #     dudx_tmp .= 0
    #     for e in 1:length(x)
    #         _, dρe = get_ρ_dρ(x[e], penalty, xmin)
    #         celldofs!(global_dofs, dh, e)
    #         Keu = bcmatrix(Kes[e]) * u[global_dofs]
    #         dudx_tmp[e] = -dρe * dot(Keu, solver.lhs[global_dofs])
    #     end
    #     return nothing, dudx_tmp # J1' * v, J2' * v