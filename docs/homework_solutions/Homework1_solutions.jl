################################################################################
# Problem 1
################################################################################

function factorial2(n::Int64)
    val = 1 #initialize value
    if n < 0
        print("You can't do that!")
        return #this returns an empty value
    else
        for i = 1:n
            val*=i #multiply by all integers between 1 and n
        end
        return val #return
    end
    
end

factorial2(3)
factorial2(-1)
#note: if n=0, setting the initial value =1 means that the correct value (1) will be returned even though 1:n is not a valid range when n=0. In this case, it just skips the loop if n=0. 

################################################################################
# Problem 2
################################################################################

function p(x::Float64, coeff::Vector{Float64})
    val = 0 #initializes
    for (i, coef) in enumerate(coeff)
        val += x^(i-1) * coef
    end
    val
end
p(1.0, [2.0, 3.0, 4.0])



################################################################################
# Problem 3
################################################################################

function approx_pi(iter::Int64)
    avg = 0
    for i = 1:iter
        x, y = rand(2) #bivariate uniform draw
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2) #distance from middle of square; pythagorean theorem
        avg+=(dist<0.5)
    end
    avg/=iter #approximation of area
    avg*=4 #π = A / r^2, where r = 1/2. So here π = 4A
    avg
end


#vectorized version of the above which eliminates the need for a loop - runs much faster
function approx_pi_v(iter::Int64)
    #create random vectors with the desired number of values
    x,y = rand(iter), rand(iter)
    #if (x,y) is within a circle of radius r around the origin, it must be the case that x^2 + y^2 <= r^2. This means we can check this by squaring the elements of both x and y, finding their sum, and taking the square root of it. 
    #the code below does this using the . syntax (which broadcasts the operation to each element of the vector/array)
    dist_from_origin = sqrt.(x.^2 + y.^2)
    #this line checks whether each element of dist_from_origin is less than 1.0. It creates a vector of length iter which is equal to 1 if that entry of dist_from_origin is less than 1, and zero otherwise. Then, to get the number of successes, it takes the sum of this vector. 
    within_circle = sum(dist_from_origin.<=1.0)
    
    avg = 4*within_circle/iter

    return avg
end
using BenchmarkTools
#median runtime around 26ms on my computer
@benchmark approx= approx_pi(1000000)

#median runtime around 6ms on my computer!
@benchmark approx = approx_pi_v(1000000)



################################################################################
# Problem 4
################################################################################

using Distributions, Plots

function simulate_ols(M::Int64, N::Int64)
    a, b, c, d, σ = 0.1, 0.2, 0.5, 1.0, 0.1 #true parameter values
    #preallocate an array to store the ols parameter results
    ols_vec = zeros(M, 4)
    x_1, x_2 = rand(Normal(), 50), rand(Normal(), 50) #random draws

    #loop over simulations
    for i = 1:M
        w = rand(Normal(), N) #randomly draw values of w

        #construct X and Y
        X = hcat(x_1, x_1.^2, x_2, ones(N))
        Y = a .* x_1 .+ b.*x_1.^2 .+ c.*x_2 .+ 1.0 .+ σ .* w
        ols = inv(X' * X) * X' * Y
        ols_vec[i,:] = ols
    end
    ols_vec
end

#Simulate and plot
ols_vec = simulate_ols(200,50)
p1 = histogram(ols_vec[:,1], title="Estimates of a")
p2 = histogram(ols_vec[:,2], title="Estimates of b")
p3 = histogram(ols_vec[:,3], title="Estimates of c")
p4 = histogram(ols_vec[:,4], title="Estimates of d")
plot(p1, p2, p3, p4, layout = (2,2), legend=:none)


################################################################################
# Problem 5
################################################################################

using Distributions, Plots

function simulate_walk(n::Int64, α::Float64, σ::Float64, t_max::Int64)
    first_times = zeros(n)

    for i = 1:n #loop over simulations
        cross_yet = false #we want to keep track of whether we have crossed 0 yet - once we do, we want to be able to keep track of this. We start out with cross_yet=false and only update it if we cross
        x_now = 1
        for t = 1:t_max
            x_next = α * x_now + σ * rand(Normal())
            if x_next<=0 #first crossing
                first_times[i] = t
                cross_yet = true #keep track of whether we crossed 0 in this iteration. This will be used to check if we ever crossed
                break #we can safely stop iterating once we reach the first value less than or equal to zero. We do this using the break command. This can save on a bunch of iterations (especially for alpha=0.8)
            end
            x_now = x_next
        end
        if ! cross_yet #if we never crossed in the first 200 iterations, set the first passage time equal to t_max
            first_times[i] = t_max
        end
    end
    first_times #return
end

#0.8
temp = simulate_walk(100, 0.8, 0.2, 200)
histogram(temp)
#to save, use the savefig command 
#1.0
temp = simulate_walk(100, 1.0, 0.2, 200)
histogram(temp)

#1.2
temp = simulate_walk(100, 1.2, 0.2, 200)
histogram(temp)



################################################################################
# Problem 6
################################################################################

function Newton(f, fp, x_0::Float64, tol::Float64, maxiter::Int64)
    error, root = 100, 0
    while error>tol
        root = x_0 - f(x_0)/fp(x_0)
        error = abs(root - x_0)
        x_0 = root
    end
    root
end

f(x) = (x-1)^3
fp(x) = 3 * (x-1)^2
root = Newton(f, fp, 4.0, 1e-3, 100)

#different function: f(x) = ln(x) + 3x - 7
f(x) = 2*log(x) + 3*x - 7
fp(x) = (2/x) + 3
root = Newton(f, fp, 4.0, 1e-3, 100)



###################################################################################################
# Question 7: Optimal Investment Problem
###################################################################################################

using Parameters 

### Struct for our model paramters
@with_kw struct ModelParameters

    β::Float64 = 0.99 
    δ::Float64 = 0.025
    α::Float64 = 0.36
    k_grid::Vector{Float64}=collect(range(0.1, length = 1000, stop = 45.0))       
    N_k::Int64 = length(k_grid)
    z_grid::Vector{Float64} = [1.25; 0.2]
    N_z::Int64 = 2
    Π::Array{Float64,2} = [0.977 0.023; 0.074 0.926]
    tol::Float64 = 10^-4

end

### Struct for our Model solutions
@with_kw struct ModelSolutions
    #because the problem has a two dimensional state vector (1000 capital levels and 2 productivity states), the value and policy functions will be arrays instead of vectors
    V::Array{Float64,2} 
    kp::Array{Float64,2}

end 

#function which constructs a ModelSolutions object with pre-specified zero guesses
function build_ModelSolutions(para)

    V = zeros(Float64,para.N_k, para.N_z)
    kp = zeros(Float64,para.N_k, para.N_z)

    sols = ModelSolutions(V,kp)

    return sols

end

function build_structs()
    #create para, a struct of type ModelParameters with the default values
    para = ModelParameters()
    #create sols, the struct of type ModelSolutions with the initial guesses
    sols = build_ModelSolutions(para)

    return para, sols

end


### Bellman operator
function bellman(para, sols)
    #unpack all of the objects within para, which is of type ModelParameters
    @unpack_ModelParameters para
    #unpack all of the objects within sols, which is of type ModelSolutions
    @unpack_ModelSolutions sols

    #create the empty initial guesses
    V_next = zeros(Float64,N_k,N_z)
    kp_next = zeros(Float64,N_k,N_z)

    #using eachindex is nice because it makes sure your loops are of the correct size and contain all of the indices of the thing you're looping over
    for i_k = eachindex(k_grid), i_z = eachindex(z_grid)
        # a bad initial value to compare against
        max_util = -1e10
         k = k_grid[i_k] #extract the capital level corresponding to index i_k
         z = z_grid[i_z] #extract the productivity level corresponding to index i_z
        budget = z*k^α + (1-δ)*k #compute the budget constraint

        for i_kp = eachindex(k_grid) #iterate over values of capital tomorrow
            kp = k_grid[i_kp] #extract the level of capital tomorrow corresponding to index i_kp
            #based on the choice of kp, find the amount of consumption left over and call it c
             c = budget - kp

            if c > 0 #if consumption is positive; otherwise log(c) will be undefined
                #the utility of choosing kp today, as well as the discounted expected value of having kp tomorrow
                #note that this is just a linear combination of the two value functions, with the weights given by their respective probabilities given the current state i_z 
                 V_temp = log(c) + β*(Π[i_z,1]*V[i_kp,1] + Π[i_z,2]*V[i_kp,2]) 
                #if the value of choosing kp today is better than the previous best value, update max_util and the policy function
                if V_temp > max_util
                    max_util = V_temp
                    kp_next[i_k,i_z] = kp
                end

            end

        end
        #update the value function for state combination i_k, i_z with the highest observed utility level (i.e. the value of capital level kp_next[i_k, i_z])
         V_next[i_k,i_z] = max_util

    end
    #return the updated value and policy functions
    return V_next, kp_next

end



### Solve model
function solve_model(para, sols)    
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    #V_next = zeros(N_k,N_z)
    #kp_next = zeros(N_k,N_z)
    max_diff = tol + 10.0
    n = 0

    while max_diff > tol
        n +=1
        V_next, kp_next = bellman(para,sols)
        
        max_diff = maximum(abs.(V_next - V))
        V .= V_next
        kp .= kp_next

        @show n, max_diff

    end

end

para, sols = build_structs();

@elapsed solve_model(para,sols)
#so much faster than the example code! On my computer, this runs in 6.16 seconds

using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))


