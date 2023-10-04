using TopOpt, Test, Zygote, Test
fc′ = 1:0.1:4.0
f2e = x-> 4700*sqrt(x)
Es = [1e-5, 1.0, 4.0] # Young's moduli of 3 materials (incl. void)
Es = vcat( [1e-5],f2e.(fc′) )

f2g = x-> 0.5*sqrt(x)/5
densities = [0.0, 0.5, 1.0] # for mass calc
densities = vcat( [0.0], f2g.(fc′) )./0.2
nmats = 3
nmats = length(Es)
nu = 0.3 # Poisson's ratio
f = 1.0 # downward force

# problem definition
problem = PointLoadCantilever(
    Val{:Linear}, # order of bases functions
    (160, 40), # number of cells
    (1.0, 1.0), # cell dimensions
    1.0, # base Young's modulus
    nu, # Poisson's ratio
    f, # load
)
ncells = TopOpt.getncells(problem)

# FEA solver
solver = FEASolver(Direct, problem; xmin=0.0)

# density filter definition
filter = DensityFilter(solver; rmin=4.0)

# compliance function
comp = Compliance(solver)

# Young's modulus interpolation for compliance
penalty1 = TopOpt.PowerPenalty(3.0) # take young modulus in each material 
interp1 = MaterialInterpolation(Es, penalty1)

# density interpolation for mass constraint
penalty2 = TopOpt.PowerPenalty(1.0) #no penalty.
interp2 = MaterialInterpolation(densities, penalty2)

# objective function
obj = y -> begin # y is a decision varible,
    x = tounit(MultiMaterialVariables(y, nmats)) 
    _E = interp1(filter(x))
    return comp(_E) #take that and multiply by the volume
end

# initial decision variables as a vector
y0 = zeros(ncells * (nmats - 1))

# testing the objective function
obj(y0)
# testing the gradient
Zygote.gradient(obj, y0)

# mass constraint
constr = y -> begin 
    _rhos = interp2(MultiMaterialVariables(y, nmats)) #rho is the density
    return sum(_rhos.x) / ncells - 0.2 # elements have unit volumes, 0.4 is the target.
end

# testing the mass constraint
constr(y0)
# testing the gradient
Zygote.gradient(constr, y0)

# building the optimization problem
model = Model(obj)
addvar!(model, fill(-10.0, length(y0)), fill(10.0, length(y0)))
add_ineq_constraint!(model, constr)

# optimization settings
alg = MMA87()
options = MMAOptions(; s_init=0.1, tol=Tolerance(; kkt=1e-3))

y0 = zeros(ncells * (nmats - 1))

# solving the optimization problem
res = optimize(model, alg, y0; options)
y = res.minimizer

# testing the solution
@test constr(y) < 1e-6

x = TopOpt.tounit(reshape(y, ncells, nmats - 1))  #reshape into a matrix with rows of each elements, and columns of each material propabilities/
sum(x[:, 2:3]) / size(x, 1) # the non-void elements as a ratio
@test all(x -> isapprox(x, 1), sum(x; dims=2))


using Makie, GLMakie
using TopOpt.TrussTopOptProblems.TrussVisualization: visualize
#create mapping from 1- 6400 to grid 1600 * 40 
x_plot =  Vector{Float64}(undef, ncells)
colx = Vector{Int64}(undef, ncells)
for i in eachindex(x_plot)
    x_plot[i], colx[i] = findmax(x[i, :])
end

map = Array{Int64,2}(undef, ncells, 2)
for i in 1:160*40
    map[i,:] = [div(i, 160) + 1, mod(i, 160)]
end
cols = [:black, :red, :blue, :green]
cols = [:black, :red, :blue, :green, :yellow, :orange, :purple, :cyan, :magenta, :brown, :pink, :gray]
col_plot = cols[colx]

f1 = Figure(resolution = (800, 600))
ax1 = Axis(f1[1, 1])
scatter!(ax1, map[:,2],map[:,1], color = colx)
f2, ax ,hm = heatmap(map[:,2],map[:,1], colx)
Colorbar(f1[1, 2])
Colorbar(f2[:,end+1], hm)
f2
f1