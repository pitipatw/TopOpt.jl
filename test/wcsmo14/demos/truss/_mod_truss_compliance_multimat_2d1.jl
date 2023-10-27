# module TrussComplianceDemo2D1
using Revise
using Zygote
using TopOpt, LinearAlgebra, StatsFuns
using Makie, GLMakie
using TopOpt.TrussTopOptProblems.TrussVisualization: visualize

begin
    Etimber = 11
    Esteel = 200
    ΔE = Esteel - Etimber 

    ρtimber = 570.0
    ρsteel = 7870.0
    Δρ = ρsteel - ρtimber 

    ρECCtimber = 0.42
    ρECCsteel = 1.45
    ΔρECC = ρECCsteel - ρECCtimber

    tmin = -8.6
    smin = 0.
    Δσmin = smin - tmin

    tmax = 0.
    smax = 345.
    Δσmax = smax - tmax
    
    η = 3

    Ee(x::Float64) = x^η*ΔE + Etimber
    ρe(x::Float64) = x^η*Δρ + ρtimber
    ρECCe(x::Float64) = x^η*ΔρECC + ρECCtimber
    σmin(x::Float64) = x^η*Δσmin + tmin
    σmax(x::Float64) = x^η*Δσmax + tmax
    
end


# 2D
node_points, elements, mats, crosssecs, fixities, load_cases = load_truss_json(
    joinpath(@__DIR__, "modtim_2d.json")
)
ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
loads = load_cases["0"]
problem = TrussProblem(
    Val{:Linear}, node_points, elements, loads, fixities, mats, crosssecs
)

xmin = 0.0001 # minimum density
x0 = fill(0.6, 2*ncells) # initial design
p = 4.0 # penalty
V = 0.4 # maximum volume fraction

solver = FEASolver(Direct, problem; xmin=xmin)
comp = TopOpt.Compliance(solver)

ts = TopOpt.TrussStress(solver)
ts(PseudoDensities(x0[1:ncells]))

function truss_stress(x)
    ae = x[1:ncells]
    xe = x[ncells+1:end]
    ee = Ee.(xe)
    σ = ts(ae.*ee)
    return σ
end

le = Vector{Float64}(undef, length(elements))
for i in eachindex(elements)
    p1 = node_points[elements[i][1]]
    p2 = node_points[elements[i][2]]
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]
    le[i] = norm([x1-x2, y1-y2])
end

function constr1(x)
    # volume fraction constraint
    ae = x[1:ncells]
    xe = x[ncells+1:end]
    return sum(ae.*le) / sum(le) - V
end

function constr2_min(x)
    @show size(x)
    ae = x[1:ncells]
    xe = x[ncells+1:end]
    minσ = σmax.(xe)
    σ = ts(PseudoDensities(ae))
    return σ .- minσ
end

function constr2_max(x)
    ae = x[1:ncells]
    xe = x[ncells+1:end] 
    maxσ = σmin.(xe)   
    σ = ts(PseudoDensities(ae))
    return σ.-maxσ
end

function obj(x)
    ae = x[1:ncells]
    xe = x[ncells+1:end]
    ρecc = ρECCe.(xe)
    return sum( ae.*le.*ρecc) 
end



x0 = rand(22)
constr2_min(x0)
g1 = Zygote.gradient(x -> sum(constr2_min(x)),x0)[1]
using FiniteDifferences
g2 = FiniteDifferences.grad(central_fdm(5, 1), x -> sum(constr2_min(x)), x0)[1]
norm(g1-g2)



m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr1)
Nonconvex.add_ineq_constraint!(m, constr2_min)
Nonconvex.add_ineq_constraint!(m, constr2_max)

options = MMAOptions(; maxiter=1000, tol=Tolerance(; kkt=1e-4, f=1e-4))
TopOpt.setpenalty!(solver, p)
@time r = Nonconvex.optimize(
    m, MMA87(; dualoptimizer=ConjugateGradient()), x0; options=options
)

@show obj(r.minimizer)
@show constr1(r.minimizer)
@show constr2_max(r.minimizer)
@show constr2_min(r.minimizer)
fig = visualize(
   problem; solver.u, topology = r.minimizer,
   default_exagg_scale=0.0
)
Makie.display(fig)


