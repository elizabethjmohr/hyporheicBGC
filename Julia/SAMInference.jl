using Turing, Distributions, QuadGK, DifferentialEquations, MCMCChains, Plots, StatsPlots
using Cubature
#using Distributed
#addprocs(4)
# @everywhere using Turing, QuadGK, DifferentialEquations
# using Zygote
# using DiffEqSensitivity
# Turing.setadbackend(:zygote)

function E(τ, α, τ₋, τ₊)
    ((τ^(-α) - τ₊^(-α))*(-α+1.0)/(τ₊^(-α+1.0)-τ₋^(-α+1.0)-τ₊^(-α)*(τ₊ - τ₋)*(1.0-α)))
end

# Specify the differential equation
function sam_model!(du,u,h,p,t)
  τ₋, τ₊, α, k, qPerUnitHZ, sHZ, sCH = p
  if(t < τ₋)
      integral = h(p, t-τ₋) * quadgk(τ -> E(τ, α, τ₋, τ₊)*exp(-k*τ),τ₋, τ₊)[1]
  elseif(τ₋<t<τ₊)
      integral = h(p, t-τ₊) * quadgk(τ -> E(τ, α, τ₋, τ₊)*exp(-k*τ), t, τ₊)[1] +
      quadgk(τ -> E(τ, α, τ₋, τ₊)* h(p,t-τ)[1]*exp(-k*τ),τ₋, t)[1]
  else
      integral = quadgk(τ -> E(τ, α, τ₋, τ₊)* h(p,t-τ)[1]*exp(-k*τ),τ₋, τ₊)[1]
  end
  du[1] = qPerUnitHZ*(sHZ/sCH)*(integral - u[1])
end

# Generate synthetic data with reasonable parameter values
τ₋ = 0.1
τ₊ = 1000.0
α = 1.6
k = 0.0
s = 27.0
sHZ = s/3
sCH = s - sHZ
qPerUnitHZ = (τ₊^(-α+1) - τ₋^(-α + 1) - (τ₊^(-α))*(τ₊-τ₋)*(-α+1))/
    ((1/2)*(τ₊^(-α))*(1-α)*(τ₊^2-τ₋^2)-(1/(2-α))*(τ₊^(2-α)-τ₋^(2-α))-(τ₊^(1-α))*(τ₊-τ₋)*(-α))
p = [τ₋, τ₊, α, k, qPerUnitHZ, sHZ, sCH]
tspan = (0.0,1000.0)
h(p,t) = 310.0
u0 = [325.0]
prob = DDEProblem(sam_model!, u0, h, tspan, p)
@time sol = solve(prob, MethodOfSteps(RK4()), saveat = 10)
fakeData = Array(sol) + 0.1 * randn(size(Array(sol)))
plot(sol, alpha = 0.3, legend = false); scatter!(sol.t, fakeData')

@model function fit_sam_model(data, prob)
    # σ ~ truncated(Cauchy(0.01,0.2),0.01,Inf)
    σ ~ Uniform(0,0.5)
    α ~ Normal(1.5,0.1)
    τ₋ = 0.01
    τ₊ = 1.0e3
    # τ₊ ~ Uniform(10,3600) #maximum RT ranging from 10s to 1 hour
    qPerUnitHZ = (τ₊^(-α+1) - τ₋^(-α + 1) - (τ₊^(-α))*(τ₊-τ₋)*(-α+1))/
    ((1/2)*(τ₊^(-α))*(1-α)*(τ₊^2-τ₋^2)-(1/(2-α))*(τ₊^(2-α)-τ₋^(2-α))-(τ₊^(1-α))*(τ₊-τ₋)*(-α))
    #α ~ truncated(Normal(1.6,0.1),1.2,2.0)
    p = [0.01, τ₊, α, 0.0, qPerUnitHZ, 9.0,18.0]
    prob = remake(prob, p = p)
    predicted = solve(prob, MethodOfSteps(RK4()), saveat= 10)
    for i = 1:length(predicted)
        data[:,i][1] ~ Normal(predicted[i][1], σ)
    end
end

model = fit_sam_model(fakeData, prob)
@time out = sample(model, MH(), 10)
plot(out)
α_summary = out[:α]
plot(α_summary, seriestype = :histogram)

# using Optim
# mle_estimate = optimize(model, MLE())

# This next command runs 3 independent chains without using multithreading.
chain = mapreduce(c -> sample(model, MH(),1000), chainscat, 1:3)

### Data import
using CSV, DataFrames
flume1Data = CSV.read("/Users/elizabethmohr/Documents/MSU/RProjects/SAM_Data_Analysis/cleanedData/SCFlume1.csv", DataFrame)

### Playing with numerical integration
using ApproxFun
integrand = function(τ)
    E(τ, 1.6, 0.1, 1000)
end
d = Interval(τ₋, τ₊)
Σ = DefiniteIntegral(Chebyshev(Interval(1,100)))
f = Fun(integrand, d)
Σ*f
