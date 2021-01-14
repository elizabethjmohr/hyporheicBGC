using DifferentialEquations
using QuadGK
using Plots

# Define model function
function sam_model!(du,u,h,p,t)
  τ₋, τ₊, α, k, q = p
  if(τ₋<t<τ₊)
    integral, error = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1)/(τ₊^(-α+1)-τ₋^(-α+1)-τ₊^(-α)*(τ₊ - τ₋)*(1-α))) * h(p,t-τ)[1]*exp(-k*τ),τ₋, t, τ₊)
  else
    integral, error = quadgk(τ -> ((τ^(-α) - τ₊^(-α))*(-α+1)/(τ₊^(-α+1)-τ₋^(-α+1)-τ₊^(-α)*(τ₊- τ₋)*(1-α))) * h(p,t-τ)[1]*exp(-k*τ), τ₋, τ₊)
  end
  du[1] = q*(integral - u[1])
end

function simulateSAM(τ₋ = 1e-3,τ₊ = 1e3, α = 1.6, k = 0.0, S = 1.0, hist = 0.0, C₀ = [4.0], t₀ = 0.0, ⏰ = 1e3)
  # Calculate hyporheic exchange (volume per time)
  q = S/quadgk(τ -> (τ₊^(-α+1.00) - τ^(-α+1.00)- (τ₊^(-α)*(τ₊-τ)*(-α+1.00)))/(τ₊^(-α+1.00) - τ₋^(-α+1.00)- (τ₊^(-α)*(τ₊-τ₋)*(-α+1.00))), τ₋, τ₊)[1]
  # Put model parameters into a vector
  p = [τ₋, τ₊, α, k, q]
  # Make a tuple containing the start and end times of the simulation
  tspan = (t₀, ⏰)
  # Specify condition before simulation start
  h(p,t) = hist
  # Build the DDE problem
  prob = DDEProblem(sam_model!, C₀, h, tspan, p)
  # Specify algorithm
  alg = MethodOfSteps(RK4())
  # Solve problem
  sol = solve(prob, alg, reltol = 1e-6)
  return(sol)
end

# Solve with default parameters
solution = simulateSAM()

# Plot solution
plot(solution)
