using DiffEqOperators, DifferentialEquations, Plots

# Automated Finite-Difference for Symbolic PDES Method
using ModelingToolkit
@parameters t x v
@variables(u(..))
@derivatives Dt'~t
@derivatives Dx'~x
v = 1.0

#1D PDE and boundary conditions
eq = Dt(u(t,x)) ~ -v*Dx(u(t,x))
bcs = [u(0,x) ~ 0,
        u(t,0) ~ 1.0,
        Dx(u(t,10.0)) ~ 0.0]

# Space and time domains
domains = [t ∈ IntervalDomain(0.0,10.0),
        x ∈ IntervalDomain(0.0,10.0)]

# PDE system
pdesys = PDESystem(eq,bcs,domains,[t,x],[u])

# Method of lines discretization
dx = 0.01
order = 3
discretization = MOLFiniteDifference(dx,order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
sol = solve(prob,Tsit5(),saveat=1.0)

# Plot results and compare with exact solution
t = sol.t
x = prob.space[2]

plt = plot()
for i in 1:length(t)
    plot!(x,Array(prob.extrapolation[1](t[i])*sol.u[i]),label="t=$(t[i])")
end
display(plt)

# Finite Difference Operator Method
tspan = (0.0, 10.0)
xₙ = 10.0 # Physical length of simulation
Δx = 0.01 # Grid spacing along x dimension (longitudinal)
N = convert(Int64, xₙ/Δx - 1) # Number of interior points
xInt = range(Δx, step = Δx,length = N)
x = range(0, step = Δx, length = (N+2))
order = 3
Δ = UpwindDifference(1,1,Δx,N,1.) # First derivative operator with respect to x
bc = RobinBC((1.0, 0.0, 1.0),(0.0, 1.0, 0.0), Δx, 1) # Boundary conditions
#bc = DirichletBC(1.0,0.0)
u₀ = zeros(N)
for i in 1:3
  u₀[i] = 1
end
p = [v]
function advect!(du,u,p,t)
        du = -p[1]*(Δ*bc*u)
end
prob = ODEProblem(advect!, u₀, tspan, p)
alg = Rosenbrock23(autodiff=false)
maxdt = Δx/v
sol = solve(prob, alg, saveat = 1.0, dtmax = 0.8*maxdt)
t = sol.t
u = sol.u

plt = plot()
for i in 1:length(t)
    plot!(x, bc*u[[i]][1] ,label="t=$(t[i])")
end
display(plt)
