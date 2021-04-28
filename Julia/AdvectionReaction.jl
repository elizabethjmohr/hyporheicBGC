using DiffEqOperators, DifferentialEquations, QuadGK
# Finite Difference Operator Method
τ₋ = 1.0 # minimum residence time
τ₊ = 100.0 # maximum residence time
α = 1.6 # power-law exponent
k = 0.0 # First-order rate constant
S = 1.0 # Volume of water in hyporheic zone
# Calculate rate of hyporheic exchange (volume per time)
q = S/quadgk(τ -> (τ₊^(-α+1.00) - τ^(-α+1.00)- (τ₊^(-α)*(τ₊-τ)*(-α+1.00)))/(τ₊^(-α+1.00) - τ₋^(-α+1.00)- (τ₊^(-α)*(τ₊-τ₋)*(-α+1.00))), τ₋, τ₊)[1]
v = 1 # Velocity
p = [τ₋, τ₊, α, k, q, v] #Vector of parameters

t₀ = 0.0 # Start of simulation time
⏰ = 5.0 # End of simulation time
tspan = (t₀, ⏰)

xₙ = 10 # Physical length of simulation
N = 999 # Number of interior points
Δx = xₙ/(N+1) # Grid spacing along x dimension (longitudinal)
x = range(0, step = Δx, length = (N+2)) # x-points of solution, including boundaries
h(p,t) = zeros(N) # History function
# Initial condition
C₀ = zeros(N)
for i in 1:N
    C₀[i] = 1.0 - i*Δx/xₙ
end
order = 4
ΔUpwind = UpwindDifference(1,order,Δx,xₙ) # First Derivative Operator, third-order approximation
ΔCentered = CenteredDifference(1,order,Δx,xₙ,1)
Robinbc = RobinBC((1.0, 0.0, 1.0),(0.0, 1.0, -0.1), Δx, order) # Boundary conditions
Dirichletbc = DirichletBC(1.0,0.0)

function advectionReaction!(du,u,p,t)
  τ₋, τ₊, α, k, q, v = p
   integral = zeros(length(u))
  if(τ₋<t<τ₊)
      integral = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1)/(τ₊^(-α+1)-τ₋^(-α+1)-τ₊^(-α)*(τ₊ - τ₋)*(1-α))) * h(p,t-τ) * exp(-k*τ),τ₋, t, τ₊)[1]
  else
      integral = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1)/(τ₊^(-α+1)-τ₋^(-α+1)-τ₊^(-α)*(τ₊- τ₋)*(1-α))) * h(p,t-τ) * exp(-k*τ), τ₋, τ₊)[1]
  end
  du = -v*ΔCentered*Dirichletbc*u #+ q*(integral - u)
end

# Build the DDE problem
prob = ODEProblem(advectionReaction!, C₀, tspan, p)
# Specify algorithm
alg = Rosenbrock23()
# Solve problem
solution = solve(prob,alg)

# Plot results and compare with exact solution
t = solution.t
u = solution.u

using Plots
plt = plot()
for i in 1:length(t)
    plot!(x, Dirichletbc*u[[i]][1] ,label="t=$(t[i])")
end
display(plt)

# Automated Finite-Difference for Symbolic PDES Method
using ModelingToolkit
@parameters t x τ₋ τ₊ α k q
@variables(u(..))
@derivatives Dt'~t
@derivatives Dx'~x

#1D PDE and boundary conditions
eq = Dt(u(t,x)) ~ -Dx(u(t,x))
bcs = [u(0,x) ~ 0,
        u(t,0) ~ 1.0,
        Dx(u(t,10.0)) ~ 0.0]

# Space and time domains
domains = [t ∈ IntervalDomain(0.0,10.0),
        x ∈ IntervalDomain(0.0,10.0)]

# PDE system
pdesys = PDESystem(eq,bcs,domains,[t,x],[C])

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
