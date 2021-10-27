using Test
using LinearAlgebra
using Ferrite
using NonconvexIpopt

using TopOpt
using Arpack

fea_ins_dir = joinpath(@__DIR__, "instances", "fea_examples");
gm_ins_dir = joinpath(@__DIR__, "instances", "ground_meshes");

# @testset "Tim buckling log-barrier $problem_dim" for problem_dim in ["2d"] # , "3d"
#     file_name = "tim_$(problem_dim).json"
#     problem_file = joinpath(gm_ins_dir, file_name)

#     mats = TrussFEAMaterial(1.0, 0.3);
#     crossecs = TrussFEACrossSec(800.0);

#     node_points, elements, _, _ , fixities, load_cases = load_truss_json(problem_file)
#     ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
#     loads = load_cases["0"]

#     problem = TrussProblem(Val{:Linear}, node_points, elements, loads, fixities, mats, crossecs);

#     xmin = 0.0001 # minimum density
#     p = 4.0 # penalty
#     V = 0.5 # maximum volume fraction
#     x0 = fill(1.0, ncells) # initial design

#     solver = FEASolver(Direct, problem);

#     # # * Before optimization, check initial design stability
#     # solver.vars = x0
#     # solver()
#     # K, G = buckling(problem, solver.globalinfo, solver.elementinfo; u=solver.u);
#     # @test isfinite(logdet(cholesky(K+G)))
#     # sparse_eigvals, buckmodes = eigs(-G,K, nev=1, which=:LR)
#     # smallest_pos_eigval = 1/sparse_eigvals[1]
#     # @test smallest_pos_eigval >= 1.0

#     Nonconvex.NonconvexCore.show_residuals[] = true

#     comp = TopOpt.Compliance(problem, solver)
#     # TODO "manual" interior point loop, adjusting the c value every iter
#     for c in [0.1] # 10:-0.1:0.1
#         function obj(x)
#             solver.vars = x
#             # trigger assembly
#             solver()
#             K, G = buckling(problem, solver.globalinfo, solver.elementinfo);
#             # minimize compliance
#             return comp(x) - c*logdet(cholesky(Array(K+G)))
#         end
#         function constr(x)
#             # volume fraction constraint
#             return sum(x) / length(x) - V
#         end

#         m = Model(obj)
#         addvar!(m, zeros(length(x0)), ones(length(x0)))
#         Nonconvex.add_ineq_constraint!(m, constr)

#         options = MMAOptions(
#             maxiter=1000, tol = Tolerance(kkt = 1e-4, f = 1e-4),
#         )
#         TopOpt.setpenalty!(solver, p)
#         r = Nonconvex.optimize(
#             m, MMA87(dualoptimizer = ConjugateGradient()),
#             x0, options = options,
#         );
#     end

#     # check result stability
#     solver = FEASolver(Direct, problem; xmin=xmin);
#     solver.vars = r.minimizer;
#     solver()
#     K, G = buckling(problem, solver.globalinfo, solver.elementinfo; u=solver.u);
#     @test isfinite(logdet(cholesky(K+G)))
#     sparse_eigvals, buckmodes = eigs(-G,K, nev=1, which=:LR)
#     smallest_pos_eigval = 1/sparse_eigvals[1]
#     @test smallest_pos_eigval >= 1.0

#     # using Makie
#     # using TopOpt.TrussTopOptProblems.TrussVisualization: visualize
#     # fig = visualize(problem; topology=r.minimizer)
# end

@testset "Tim buckling SDP constraint $problem_dim" for problem_dim in ["2d"] # , "3d"
    file_name = "tim_$(problem_dim).json"
    problem_file = joinpath(gm_ins_dir, file_name)

    mats = TrussFEAMaterial(1.0, 0.3);
    crossecs = TrussFEACrossSec(800.0);

    node_points, elements, _, _ , fixities, load_cases = load_truss_json(problem_file)
    ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
    loads = load_cases["0"]

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, fixities, mats, crossecs);

    xmin = 0.0001 # minimum density
    p = 4.0 # penalty
    V = 0.5 # maximum volume fraction
    x0 = fill(V, ncells) # initial design

    solver = FEASolver(Direct, problem);
    ch = problem.ch
    dh = problem.ch.dh

    comp = TopOpt.Compliance(problem, solver)
    dp = TopOpt.Displacement(solver)
    assemble_k = TopOpt.AssembleK(problem)
    element_k = ElementK(solver)
    truss_element_kσ = TrussElementKσ(problem, solver)

    # * comliance minimization objective
    obj = comp
    c = 1.0 # buckling load multiplier

    function buckling_matrix_constr(x)
        # * Array(K + c*Kσ) ⋟ 0, PSD
        # * solve for the displacement
        u = dp(x)

        # * x -> Kes, construct all the element stiffness matrices
        # a list of small matrices for each element (cell)
        Kes = element_k(x)

        # * Kes -> K (global linear stiffness matrix)
        K = assemble_k(Kes)
        K = apply_boundary_with_meandiag!(K, ch)

        # * u_e, x_e -> Ksigma_e
        Kσs = truss_element_kσ(u, x)

        # * Kσs -> Kσ
        Kσ = assemble_k(Kσs)
        Kσ = apply_boundary_with_zerodiag!(Kσ, ch)

        return Array(K + c*Kσ)
    end

    function vol_constr(x)
        # volume fraction constraint
        return sum(x) / length(x) - V
    end

    # * Before optimization, check initial design stability
    @test isfinite(logdet(cholesky(buckling_matrix_constr(x0))))

    m = Model(obj)
    addvar!(m, zeros(length(x0)), ones(length(x0)))
    Nonconvex.add_ineq_constraint!(m, vol_constr)
    Nonconvex.add_sd_constraint!(m, buckling_matrix_constr)

    Nonconvex.NonconvexCore.show_residuals[] = true
    alg = SDPBarrierAlg(sub_alg=IpoptAlg())
    options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
    r = Nonconvex.optimize(m, alg, x0, options = options)
    # println("$(r.convstate)")

    # * check result stability
    @test isfinite(logdet(cholesky(buckling_matrix_constr(r.minimizer))))

    # using Makie
    # using TopOpt.TrussTopOptProblems.TrussVisualization: visualize
    # fig = visualize(problem; topology=r.minimizer)
end