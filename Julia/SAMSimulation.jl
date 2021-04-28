using DifferentialEquations
using QuadGK
using Plots

# Define model function
function sam_model!(du,u,h,p,t)
  τ₋, τ₊, α, k, q = p
  if(τ₋<t<τ₊)
    integral, error = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1.0)/(τ₊^(-α+1.0)-τ₋^(-α+1.0)-τ₊^(-α)*(τ₊ - τ₋)*(1.0-α))) * h(p,t-τ)[1]*exp(-k*τ),τ₋, t, τ₊)
  else
    integral, error = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1.0)/(τ₊^(-α+1.0)-τ₋^(-α+1.0)-τ₊^(-α)*(τ₊- τ₋)*(1.0-α))) * h(p,t-τ)[1]*exp(-k*τ), τ₋, τ₊)
  end
  du[1] = q*(integral - u[1])
end

τ₋ = 1.0e-3
τ₊ = 1.0e3
α = 1.6
k = 0.0
S = 1.0
hist = 0.0
u0 = [4.0]
t₀ = 0.0
⏰ = 1.0e3
q = S/quadgk(τ -> (τ₊^(-α+1.00) - τ^(-α+1.00)- (τ₊^(-α)*(τ₊-τ)*(-α+1.00)))/(τ₊^(-α+1.00) - τ₋^(-α+1.00)- (τ₊^(-α)*(τ₊-τ₋)*(-α+1.00))), τ₋, τ₊)[1]
p = [τ₋, τ₊, α, k, q]
tspan = (t₀, ⏰)
h(p,t) = zeros(1)
prob = DDEProblem(sam_model!, u0, h, tspan, p)
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg, saveat = 10)
fakeData = Array(sol) + 0.1 * randn(size(Array(sol)))
plot(sol, alpha = 0.3, legend = false); scatter!(sol.t, fakeData')

function simulateSAM(;τ₋ = 1e-3,τ₊ = 1e3, α = 1.6, k = 0.0, S = 1.0, hist = 0.0, C₀ = [4.0], t₀ = 0.0, ⏰ = 1e3)
  # Calculate hyporheic exchange (volume per time)
  q = S/quadgk(τ -> (τ₊^(-α+1.00) - τ^(-α+1.00)- (τ₊^(-α)*(τ₊-τ)*(-α+1.00)))/(τ₊^(-α+1.00) - τ₋^(-α+1.00)- (τ₊^(-α)*(τ₊-τ₋)*(-α+1.00))), τ₋, τ₊)[1]
  # Put model parameters into a vector
  p = [τ₋, τ₊, α, k, q]
  # Make a tuple containing the start and end times of the simulation
  tspan = (t₀, ⏰)
  # Specify condition before simulation start
  h(p,t) = hist
  # Build the DDE problem
  prob = DDEProblem(sam_model, C₀, h, tspan, p)
  # Specify algorithm
  alg = MethodOfSteps(RK4())
  # Solve problem
  sol = solve(prob, alg)
  return(sol)
end


# Solve with default parameters
lowAlpha = simulateSAM(k = 0.05, τ₋ = 1e-3,τ₊ = 1e3, α = 1.3, ⏰ = 1e3)
highAlpha = simulateSAM(k = 0.005, τ₋ = 1.00e-3,τ₊ = 1.00e3, α = 1.9, ⏰ = 1.00e3)
zeroK = simulateSAM(τ₋ = 1.00e-3,τ₊ = 1.00e3, α = 1.3, ⏰ = 1.00e3)
zeroK2 = simulateSAM(τ₋ = 1.00e-3,τ₊ = 1.00e3, α = 1.9, ⏰ = 1.00e3)
# Plot solution
plot(lowAlpha)
plot!(highAlpha)
plot!(zeroK)
plot!(zeroK2)
