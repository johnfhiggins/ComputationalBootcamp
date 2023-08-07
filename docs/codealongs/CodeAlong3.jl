################################################################################
# Load Packages
################################################################################

using Optim, Plots

################################################################################
# Univariate Boxed Constrained Optimization
################################################################################

f(x, y) = (x-y)^2


####Brent's method is the default for univariate boxed-constrained optimization
opt = optimize(x->f(x, 1.0),   #function
                -5.0,           #lower bound
                5.0             # upper bound
                )

opt.minimizer
opt.minimum

#The "x->f(x)" operator defines a function inline that maps x to f(x)
# this is useful for when optimizing over a subset of parameters. 
# Here, we optimize over x for a fixed y. 

################################################################################
# Multivariate Optimization Optimization Derivative Free
################################################################################

#####Rosenbrock
function Rosenbrock(x::Vector{Float64})
    val = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
    val
end

##evaluate function at a bunch of points so we can plot it
x_grid = collect(-3:0.01:3)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:nx
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Rosenbrock(guess)
end

surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50,50))
contourf(x_grid, x_grid, log.(1 .+ z_grid), 
            seriescolor=:inferno, xlabel="x_1", ylabel = "x_2")
scatter!([1],[1], color = "white", label="Minimizer")

#Nelder-Mead is default option for multivariate optimization
guess = [0.0, 0.0]
opt = optimize(Rosenbrock, guess)
opt.minimizer #Success!
opt.minimum

###Let's specify LBFGS to be the solver.
opt = optimize(Rosenbrock, [-5.0;-5.0], LBFGS())
opt.minimizer #also worked
opt.minimum



################################################################################
# Multivariate Optimization Optimization with Derivates
################################################################################

#Define Gradient of the Rosenbrock function
function g(G, x::Vector{Float64})
    G[1] = -2.0 * (1.0 - x[1]) -400.0 * x[1] * (x[2] - x[1]^2)
    G[2] = 200.0 * (x[2] - x[1]^2)
    G #return
end

#Rosenbrock's Hessian
function h(H, x::Vector{Float64})
    H[1,1] = 2 - 400.0*x[2] + 1200.0*x[1]^2
    H[1,2] = -400.0 * x[1]
    H[2,1] = -400.0 * x[2]
    H[2,2] = 200.0
    H #return
end

#Newton's Method is default when providing gradient and Hessian
guess = [0.0, 0.0]
opt = optimize(Rosenbrock, g, h, guess)
opt.minimizer #works!
opt.minimum

################################################################################
# Many local minimum
################################################################################

function Greiwank(x::Array{Float64,1})
    val = (1/4000)*sum(x.^2) - prod(cos.(x./sqrt(length(x)))) + 1
    val
end

##evaluate function at a bunch of points
x_grid = collect(-5:0.01:5)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:nx
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Greiwank(guess)
end

##plots
surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50,70))
contourf(x_grid, x_grid, z_grid, seriescolor=:inferno)

#global optimum at (0,0)
guess_init = [3.0, 3.0]
opt = optimize(Greiwank, guess_init) #this fails!
opt.minimizer
opt.minimum 

#now this works!
guess_init = [2.0, 2.0]
opt = optimize(Greiwank, guess_init) #this works!
opt.minimum
opt.minimizer

#try multiple starts!
function Multistart()
    x_grid = collect(-5:2.0:5)
    nx = length(x_grid)
    minimum, minimizers = 100, [100, 100] #preallocate bad values for minimum and minimizes

    for i = 1:nx, j = 1:nx
        guess = [x_grid[i], x_grid[j]] #starting guess
        opt = optimize(Greiwank, guess) #nelder-mead with new starting guess
        if opt.minimum<minimum #new minimum!
            minimum = opt.minimum #update
            minimizers = opt.minimizer #update
        end
    end
    minimum, minimizers #return
end
min, minimizers = Multistart()



################################################################################
# OLS Example
################################################################################

using Distributions, Random

##run the same OLS as first class
dist = Normal(0,1)
β_0 = 1.0
β_1 = 2.0
β_2 = 3.0
n = 10000
x = rand(n).*10
x2 = x.^2
Random.seed!(1234)  ###important to remember!!
ϵ = rand(dist, n)
Y_true = β_0 .+ β_1.*x + β_2.*x2 .+ ϵ
X = hcat(ones(n), x, x2)
β_ols = inv(X' * X) * X' * Y_true

####Bad Squared function: This will draw new ϵ every iteration
function sq_error(β::Array{Float64,1})
    β_0, β_1, β_2 = β[1], β[2], β[3] #unpack β
    Random.seed!(1234) #Must uncomment this to get it to work
    ϵ = rand(dist, n) #draw epsilons
    Y_true = 1.0 .+ 2.0.*x + 3.0.*x2 .+ ϵ
    Y_predict = β_0 .+ β_1.*x + β_2.*x2 #true and predicted value
    error = sum((Y_true.-Y_predict).^2) #sum of squared error
    error #return
end


#do OLS with nelder-mead
guess_init = [0.0, 0.0, 0.0]
opt = optimize(sq_error, guess_init) #it works

opt.minimizer