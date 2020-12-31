using DifferentialEquations
using QuadGK
using Plots

# Specify parameters
tauMin = 1.00e-3
tauMax = 1.00e3
a = 1.6
k = 0.00
S = 1.00
q = S/quadgk(tau -> (tauMax^(-a+1.00) - tau^(-a+1.00)- (tauMax^(-a)*(tauMax-tau)*(-a+1.00)))/(tauMax^(-a+1.00) - tauMin^(-a+1.00)- (tauMax^(-a)*(tauMax-tauMin)*(-a+1.00))), tauMin, tauMax)[1]
p = [tauMin, tauMax, a, k, q]

# Specify initial condition
u0 = [4.000]

# Specify condition before simulation start
h(p,t) = 0.000

# Define model function
function sam_model!(du,u,h,p,t)
  tauMin, tauMax, a, k, q = p
  if(tauMin<t<tauMax)
    integral, error = quadgk(tau -> ((tau^(-a) - tauMax^(-a))*(-a+1)/(tauMax^(-a+1)-tauMin^(-a+1)-tauMax^(-a)*(tauMax - tauMin)*(1-a))) * h(p,t-tau)[1]*exp(-k*tau), tauMin, t, tauMax)
  else
    integral, error = quadgk(tau -> ((tau^(-a) - tauMax^(-a))*(-a+1)/(tauMax^(-a+1)-tauMin^(-a+1)-tauMax^(-a)*(tauMax - tauMin)*(1-a))) * h(p,t-tau)[1]*exp(-k*tau), tauMin, tauMax)
  end
  du[1] = q*(integral - u[1])
end
# Specify simulation time
tspan = (0.00, 1.00e3)

# Build the DDE problem
prob = DDEProblem(sam_model!, u0, h, tspan, p)

# Specify algorithm
alg = MethodOfSteps(RK4())

# Solve problem
sol = solve(prob, alg)

# Plot solution
plot(sol)
