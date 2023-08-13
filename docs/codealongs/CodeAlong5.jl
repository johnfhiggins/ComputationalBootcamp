################################################################################
# Packages
################################################################################

using Random, Plots, Distributions, Statistics, Parameters


################################################################################
# Static Arrays
################################################################################

# In the homework 1 question 7, you had to calculate the expectation w.r.t the 
# transition matrix Π. 

# You could do this via matrix multiplication: EV = Π[1,:]' * V. Here, V = [V_g; V_b].

# Working with arrays such as matrices can be quite slow because although the number of
# dimensions is fixed, the size of each dimension can change. 

using BenchmarkTools

function add_up(M,R) # M is a 2x2 matrix, R is 2x1 vector 
    
    total = 0.0
    
    for i = 1:1000000
        total += M[1,:]' * R
        total += M[2,:]' * R
    end
    
    total

end

function add_up2(M,R) # M is a 2x2 matrix; R is two element vector 
    total = 0.0

    for i = 1:1000000
        total += M[1,1] * R[1] + M[1,2] * R[2]
        total += M[2,1] * R[1] + M[2,2] * R[2]
    end
    
    total

end


M = rand(2,2)
R = rand(2) 

@benchmark add_up(M,R) # Slow
@benchmark add_up2(M,R)

M_static = SMatrix{2,2}(rand(2,2)) #This is a static matrix. It must always be 2x2
R_static = SVector{2}(R) #This is a static vector. It must always have 2 elements.  

@benchmark add_up(M_static,R) # A lot faster
@benchmark add_up(M_static,R_static) # Wow!


# Static arrays are generally good for "small" arrays. This is static arrays are stored
# in a different part of memory, and it takes a longer time to generate them. 
# For example, the following line of code on my computer never ran in a reasonable time.
# M_big_static = SMatrix{100,100}(rand(100,100))
# In general, use static arrays when they are "Small". You might have to play around 
# with this. A very rough rule of thumb is that you should consider using a normal 
# Array for arrays larger than 100 elements. 



################################################################################
# Birthday problem
################################################################################

function birthday(n, sims)
    results = zeros(sims) #preallocate monte-carlo results vector
    for i = 1:sims #loop over simulations
        days = rand(1:365, n) #draw birthdays
        results[i] = length(unique(days)) #fill in vector
    end
    results #return
end
res_20= birthday(20, 100000)
histogram(res_20)

res_50 = birthday(50, 100000) 
histogram(res_50)

res_70 = birthday(70, 100000)
histogram(res_70)



################################################################################
# Average distance between two random points in a cube
################################################################################

function Point_distance(sims)
    results = zeros(sims)
    for i = 1:sims #loop over simulations
        p1, p2 = rand(3), rand(3) #two points!
        results[i] = sqrt(sum((p1.-p2).^2))
    end
    return mean(results)
end
Point_distance(10000)



################################################################################
# Expected Value of College Given Wage offer Shock
################################################################################

@with_kw struct Primitives
    β_0::Float64 = 2.7 #wage constant
    β_1::Float64 = 0.47 #college premium
    σ::Float64 = 0.597 #wage offer SD
    α::Float64 = 1.0 #leisure
    B::Float64 = 5.0 #base consumption
    d::Float64 = 0.25 #Opportunity Cost of College
end

mutable struct Results
    emax::Array{Float64,1} #first for no college, second for college
    lfp::Array{Float64,1} #lfp probabilities
    ewage::Array{Float64,1} #Average wages
    ewage_obs::Array{Float64,1} #Average observed wages
end

#initializes model primitives and executes solution
function Solve_model(sims)
    prim = Primitives() #initialize primitives
    res = Results(zeros(2), zeros(2), zeros(2), zeros(2)) #initialize resutls
    Compute_emax(prim, res, sims) #solve model
    prim, res #return
end

function Compute_emax(prim, res, sims)
    #housekeeping
    @unpack β_0, β_1, σ, α, B, d = prim
    dist = Normal(0, σ)
    val_nwork = α + log(B)
    utils, lfps, wages = zeros(2, sims), zeros(2, sims), zeros(2,sims)

    for s = 0:1 #loop over schooling levels
        for i = 1:sims #loop over simulations
            ε = rand(dist) #draw shock and compute resultant wage
            wage = exp(β_0 + β_1*s + ε)
            util = max(log(wage), val_nwork) #max utility
            utils[s+1, i] = util  #update
            lfps[s+1, i] = (log(wage)>val_nwork)
            wages[s+1, i] = wage #update working decision
        end
    end

    res.emax[1], res.emax[2] = mean(utils[1,:]), mean(utils[2,:]) #expected max utilities
    res.lfp[1], res.lfp[2] = mean(lfps[1,:]), mean(lfps[2,:]) #LFP probabilities
    res.ewage[1], res.ewage[2] = mean(wages[1,:]), mean(wages[2,:]) #Average Wages
    res.ewage_obs[1], res.ewage_obs[2] = mean(wages[1,:].*(lfps[1,:]))/res.lfp[1], mean(wages[2,:].*(lfps[2,:]))/res.lfp[2] #Average observed Wages
end

prim, res = Solve_model(100) #run the code

#Observed College Wage premium is biased down because non-college less likely to work!
res.ewage[2] / res.ewage[1] 
res.ewage_obs[2] / res.ewage_obs[1]




################################################################################
# Quadrature
################################################################################

using FastGaussQuadrature, LinearAlgebra

#Gauss-Legendre Quadrature for approximating integral between -1 and 1
#Works well with n nodes if function is well approximated by polynomial of degree n
#Get nodes x and weights w
x, w = gausslegendre(3)

#Function F
F(x) = x^4 - 3x^3 - 4x^2
#Derivative of function f
f(x) = 4x^3 - 9x^2 - 8x
#Calculate intergral appoximation
I = sum(w .* f.(x))
#Calculate true value
F(1) - F(-1)


# Gauss-Hermite Quadrature for approximating integral of function g(x) = f(x)*exp(-x^2)
# Again, works well with n nodes if f is well approximated by polynomial of degree n
# Here we will approximate the expection E[f(x)] when x is distributed according to standard normal

#Get nodes and weights
x,w = gausshermite(3)
#Change of variable: https://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
g_tilde(x,f) = f(sqrt(2)x)/sqrt(π)

#Expectation approximation
sum(w.*g_tilde.(x,f))

#Monte Carlo Simulation Approximation
val = 0
N = 1000000
for i = 1:N
    val += f(rand(Normal()))/N
end
val




################################################################################
# Gibbs Sampler for Correlated Normal
################################################################################

# Note: for illustrative purpose only. Improve this by:
# changing mean and variance beyond (0,1)
# try different seeds then calculate average

function bigibbs(T, rho)
    x = zeros(T+1)
    y = zeros(T+1)
    for t = 1:T
        x[t+1] = randn() * sqrt(1-rho^2) + rho*y[t]
        y[t+1] = randn() * sqrt(1-rho^2) + rho*x[t+1]
    end
    return x, y
end

x,y = bigibbs(100000, 0.8)
x = x[10000:end];
y = y[10000:end];
using Statistics
mean(x), var(x), mean(y), var(y), cov(x,y)