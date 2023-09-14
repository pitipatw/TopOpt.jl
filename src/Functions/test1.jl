using UnPack


@unpack σ, transf_matrices, u_fn = ts
@unpack global_dofs, solver = u_fn
@unpack penalty, problem, xmin = solver
@unpack Kes = solver.elementinfo



ts(PseudoDensities(x0))
getA(sp::TrussProblem) = [cs.A for cs in sp.truss_grid.crosssecs]

dh = problem.ch.dh
u  = u_fn(PseudoDensities(x0))
As = getA(problem)
n_dofs = length(global_dofs)






