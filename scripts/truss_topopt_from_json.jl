using GLMakie
using TopOpt
using TopOpt.TrussTopOptProblems.TrussVisualization: visualize

ins_dir = joinpath(@__DIR__, "..", "data");
file_name = "cantilever_topopt.json"
problem_file = joinpath(ins_dir, file_name)
lc_ind = 0

node_points, elements, mats, crosssecs, fixities, load_cases = load_truss_json(
    problem_file
)
loads = load_cases[string(lc_ind)]

problem = TrussProblem(
    Val{:Linear}, node_points, elements, loads, fixities, mats, crosssecs
)

ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)

V = 0.3 # volume fraction
xmin = 0.001 # minimum density
rmin = 4.0 # density filter radius

penalty = TopOpt.PowerPenalty(1.0) # 1
solver = FEASolver(Direct, problem; xmin=xmin, penalty=penalty)
## call solver to trigger assemble!
solver()

# * Compliance
comp = TopOpt.Compliance(problem, solver)
obj = comp
volfrac = TopOpt.Volume(problem, solver)
constr = x -> volfrac(x) - V

options = MMAOptions(; maxiter=3000, tol=Nonconvex.Tolerance(; kkt=0.001))
convcriteria = Nonconvex.KKTCriteria()
x0 = fill(V, length(solver.vars))
nelem = length(x0)

m = Model(obj)
addvar!(m, zeros(nelem), ones(nelem))
add_ineq_constraint!(m, constr)

TopOpt.setpenalty!(solver, penalty.p)
result = Nonconvex.optimize(m, MMA87(), x0; options=options)

println("="^10)
println("$(result.convstate)")

fig = visualize(
    problem, topology = result.topology, vector_arrowsize = 0.1,
    vector_linewidth=0.8, default_exagg_scale=ndim == 3 ? 1.0 : 0.01,
    exagg_range = ndim == 3 ? 10.0 : 0.1,
)