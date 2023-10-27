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
'u': Displacement vectors

# Returns
stress vector `σ`, compressive stress < 0, tensile stress > 0
"""
function (ts::TrussStress{T})(u::DisplacementResult) where {T}
    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    dh = getdh(problem)
    As = getA(problem)
    @unpack Kes = solver.elementinfo
    for e in 1:length(As)
        # Ke = R' * K_local * R
        # F = R * (R' * K_local * R) * u
        celldofs!(global_dofs, dh, e)
        σ[e] = -(transf_matrices[e] * Kes[e] * u.u[global_dofs])[1] / As[e]
    end
    return copy(σ)
end

function (ts::TrussStress{T})(x::PseudoDensities) where {T}
    @unpack u_fn = ts
    ts.fevals += 1
    u = u_fn(x)
    return ts(u)
end

println("HI")

function ChainRulesCore.rrule(
    ts::TrussStress{T}, u::DisplacementResult) where {T}

    @unpack σ, transf_matrices, u_fn = ts
    @unpack global_dofs, solver = u_fn
    @unpack penalty, problem, xmin = solver
    @unpack Kes = solver.elementinfo
    # K = ak(Kes)
    dh = problem.ch.dh
    
    # Forward-pass
    σ = ts(u)
    # n_dofs = length(global_dofs)
    getA(sp::TrussProblem) = [cs.A for cs in sp.truss_grid.crosssecs]
    As = getA(problem)
    function truss_stress_pullback(Δσ)
#       #gradient will be the same size as elements, or stress itself.
        # Δσ = Vector{Float64}(undef, length(σ))
        Δu = zeros(Float64,length(u))
        # dudAs = jacobian(u, x)
        for e in 1:length(Δσ)
            celldofs!(global_dofs, dh, e)
            Δu[global_dofs] += -(transf_matrices[e] * Kes[e])[1,:] ./ As[e] .* Δσ[e]
        end
    return NoTangent(), Tangent{typeof(u)}(;u = Δu)
end
    return σ , truss_stress_pullback
end
