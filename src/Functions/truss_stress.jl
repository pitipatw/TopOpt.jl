@params mutable struct TrussStress{T} <: AbstractFunction{T}
    σ::AbstractVector{T} # stress vector, axial stress per cell
    u_fn::Displacement
    transf_matrices::AbstractVector{<:AbstractMatrix{T}}
    fevals::Int
    maxfevals::Int
end

function Base.show(::IO, ::MIME{Symbol("text/plain")}, ::TrussStress)
    return println("TopOpt truss stress function")
end

"""
    TrussStress(solver; maxfevals=10^8)

Construct the TrussStress function struct.
"""
function TrussStress(solver::AbstractFEASolver; maxfevals=10^8)

    T = eltype(solver.u)
    dim = TopOptProblems.getdim(solver.problem)
    dh = solver.problem.ch.dh
    N = getncells(dh.grid)
    σ = zeros(T, N)
    transf_matrices = Matrix{T}[]
    u_fn = Displacement(solver; maxfevals)
    R = zeros(T, (2, 2 * dim))
    for (cellidx, cell) in enumerate(CellIterator(dh))
        u, v = cell.coords[1], cell.coords[2]
        # R ∈ 2 x (2*dim)
        R_coord = compute_local_axes(u, v)
        fill!(R, 0.0)
        R[1, 1:dim] = R_coord[:, 1]
        R[2, (dim + 1):(2 * dim)] = R_coord[:, 2]
        push!(transf_matrices, copy(R))
    end
    return TrussStress(σ, u_fn, transf_matrices, 0, maxfevals)
end

"""
# Arguments
`x` = design variables

# Returns
displacement vector `σ`, compressive stress < 0, tensile stress > 0
"""
function (ts::TrussStress{T})(x::PseudoDensities) where {T}
    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    dh = getdh(problem)
    ts.fevals += 1
    u = u_fn(x)
    As = getA(problem)
    @unpack Kes = solver.elementinfo
    for e in 1:length(x)
        # Ke = R' * K_local * R
        # F = R * (R' * K_local * R) * u
        celldofs!(global_dofs, dh, e)
        σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
    end
    return copy(σ)
end


function ChainRulesCore.rrule(
    ts::TrussStress{T}, x::PseudoDensities) where {T}
    #or?
    # ts::TrussStress, x::PseudoDensities)
    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    @unpack Kes = solver.elementinfo
    # K = ak(Kes)
    dh = problem.ch.dh
    
    # Forward-pass
    σ = ts(x)
    u = u_fn(x)
    # n_dofs = length(global_dofs)
    getA(sp::TrussProblem) = [cs.A for cs in sp.truss_grid.crosssecs]
    As = getA(problem)
    function truss_stress_pullback(Δ)
#         #gradient will be the same size as elements, or stress itself.
        Δσ = Vector{Float64}(undef, length(σ))
        dudAs = jacobian(u, x)
        for e in 1:length(Δσ)
            celldofs!(global_dofs, dh, e)
            #σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
            # @show global_dofs
            # @show u[global_dofs]
            # @show transf_matrices[e]
            # @show (transf_matrices[e] * Kes[e]) / As[e]
            Δσ[e] = norm((transf_matrices[e] * Kes[e]) / As[e] * duAs[global_dofs])
        end
#             _, dρe = get_ρ_dρ(x.x[e], penalty, xmin)
#             celldofs!(global_dofs, dh, e)
#             Keu = bcmatrix(Kes[e]) * u.u[global_dofs]
#             dσdx_tmp[e] = -dρe * dot(Keu, solver.lhs[global_dofs])
#         end 

    return Tangent{typeof(ts)}(;
        σ = Δσ,
        u_fn = u ,
        transf_matrices = NoTangent(),
        fevals = NoTangent(),
        maxfevals = NoTangent()), NoTangent()
end
    return σ , truss_stress_pullback
end
# TODO complete
# """
# rrule for autodiff.
    
# du/dxe = -K^-1 * dK/dxe * u
# d(u)/d(x_e) = - K^-1 * d(K)/d(x_e) * u
#             = - K^-1 * (Σ_ei d(ρ_ei)/d(x_e) * K_ei) * u
#             = - K^-1 * [d(ρ_e)/d(x_e) * K_e * u]
# d(u)/d(x_e)' * Δ = -d(ρ_e)/d(x_e) * u' * K_e * (K^-1 * Δ)
# """
# function ChainRulesCore.rrule(
#     ::typeof(ts), x::PseudoDensities)
#     #or?
#     ts::TrussStress, x::PseudoDensities)
#     @unpack σ, transf_matrices, u_fn = ts
#     @unpack global_dofs, solver = u_fn
#     @unpack penalty, problem, xmin = solver
#     @unpack Kes = solver.elementinfo
    
#     dh = getdh(problem)
    
#     # Forward-pass
#     truss_stress = ts(x)
#     n_dofs = length(global_dofs)
#     function truss_stress_pullback(Δ)
# #         #gradient will be the same size as elements, or stress itself.
#         Δσ = Vector{Float64}(undef, length(σ))
# #         for e in 1:length(σ)
# #             celldofs!(global_dofs, dh, e)
# #             #σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
# #             Δσ[e] = Tangent{typeof(ts)}(; -(transf_matrices[e] * Kes[e])[1] / As[e])

# #             _, dρe = get_ρ_dρ(x.x[e], penalty, xmin)
# #             celldofs!(global_dofs, dh, e)
# #             Keu = bcmatrix(Kes[e]) * u.u[global_dofs]
# #             dσdx_tmp[e] = -dρe * dot(Keu, solver.lhs[global_dofs])
# #         end 

#     σ::AbstractVector{T} # stress vector, axial stress per cell
#     u_fn::Displacement
#     transf_matrices::AbstractVector{<:AbstractMatrix{T}}
#     fevals::Int
#     maxfevals::Int
#             Δ = Tangent{typeof(ts)}(;
#                 σ = 
#                 u_fn = 
#                 transf_matrices =
#                 fevals = NoTangent(),
#                 maxfevals = NoTangent())
#         return Δ, NoTangent()
#     end
#     return σ , truss_stress_pullback
# end

#     return truss_stress, Δ -> begin # v
#         solver.rhs .= Δ
#         solver(reuse_chol = true, assemble_f = false)
#         dudx_tmp .= 0
#         for e in 1:length(x)
#             _, dρe = get_ρ_dρ(x[e], penalty, xmin)
#             celldofs!(global_dofs, dh, e)
#             Keu = bcmatrix(Kes[e]) * u[global_dofs]
#             dudx_tmp[e] = -dρe * dot(Keu, solver.lhs[global_dofs])
#         end
#         return nothing, dudx_tmp # J1' * v, J2' * v
#     end
# end
