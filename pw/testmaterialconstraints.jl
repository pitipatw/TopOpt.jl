using TopOpt, LinearAlgebra, StatsFuns
using TopOpt.TrussTopOptProblems.TrussVisualization: visualize
using Makie, GLMakie

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
    ρECCe(x::Float64) = x^η*ΔρECC + ρECCCtimber
    σmin(x::Float64) = x^η*Δσmin + tmin
    σmax(x::Float64) = x^η*Δσmax + tmax
    
end

    

# xe = 0:0.01:1
# σ1 = σmin.(xe)
# σ2 = σmax.(xe)
# f1 = Figure(resolution = (1080, 720))
# ax1 = f1[1,1] = Axis(f1,
#     title = "σ1 test",
# )
# ax2 = f1[1,2] = Axis(f1,
#     title ="σ2 test")
# scatter!(ax1, xe, σ1)
# scatter!(ax2, xe, σ2)
# f1

# 2D
ndim = 2
node_points, elements, mats, crosssecs, fixities, load_cases = load_truss_json(
    joinpath(@__DIR__,"tim_$(ndim)d.json")
)
ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
loads = load_cases["1"]
problem = TrussProblem(
    Val{:Linear}, node_points, elements, loads, fixities, mats, crosssecs
)

L = Vector{Float64}(undef, length(elements))
#get list of length
for e in keys(elements)
    p1 = elements[e][1]
    p2 = elements[e][2]
    (x1,y1) = node_points[p1]
    (x2,y2) = node_points[p2]
    L[e] = norm([x1 y1] - [x2 y2])
end

xmin = 0.0001 # minimum density
x0 = fill(1.0, (ncells*2)) # initial design, 2 columns for 2 Ae and Xe
p = 4.0 # penalty
# V = 0.5 # maximum volume fraction

solver = FEASolver(Direct, problem; xmin=xmin)
comp = TopOpt.Compliance(solver)

ts = TopOpt.TrussStress(solver)
function truss_stress(x)
    ae = x[1:ncells]
    xe = x[ncells+1:]
    ee = Ee.(xe)
    σ = ts(ae.*ee)
    return σ
end

# function obj(x)
#     # minimize compliance
#     Ae = x[:,1]
#     return comp(PseudoDensities(Ae))
# end
function constr1(x)
    # volume fraction constraint
    ae = x[1:ncells]
    xe = x[ncells+1:]

    return sum(ae.*L) / sum(L) - V
end

function constr2_min(x)
    Xe = x[:,2]
    min = σmax.(Xe)
    σ = truss_stress(x)
    return σ .- min
end

function constr2_max(x)
    Xe = x[:,2]
    max = σmax.(Xe)
    σ = truss_stress(x)
    return σ.-max
end

function obj(x)
    ae = x[:,1]
    xe = x[:,2]
    ρecc = ρECCe.(xe)
    return sum( ae.*le.*ρecc) 
end

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
@show constr(r.minimizer)
fig = visualize(
   problem; solver.u, topology = r.minimizer,
   default_exagg_scale=0.0
)
Makie.display(fig)

end
