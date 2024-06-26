

################################################################################
# Packages
################################################################################

using Optim, Interpolations, Plots, Polynomials, Parameters



################################################################################
# Runge's function and phenomenon
################################################################################

Runge(x) = 1 /(1 + 25*x^2)
function Runge_2(x)
    val = 1/(1+25*x^2)
    return val
end
x_fine = collect(-1.0:0.01:1.0)
x_coarse = collect(range(-1.0, length = 6, stop = 1.0))
plot(x_fine, Runge.(x_fine))
#The ! adds the plot/scatter on top of the previous plot rather than make a new plot
scatter!(x_coarse, Runge.(x_coarse))

#5-degree polynomial fit
poly_fit = fit(x_coarse, Runge.(x_coarse), 5)
plot!(x_fine, poly_fit.(x_fine))

#try a higher-degree -- surely that will help!
x_coarse_2 = collect(range(-1.0, length = 10, stop = 1.0))
poly_fit_2 = fit(x_coarse_2, Runge.(x_coarse_2), 9)
plot!(x_fine, poly_fit_2.(x_fine))

# Let's Try Interpolation Instead
x_range = -1.0:0.25:1.0
runge_range = Runge.(x_range)
runge_linear = LinearInterpolation(x_range, runge_range)
runge_spline = CubicSplineInterpolation(x_range, runge_range)
plot(x_fine, 
    [Runge.(x_fine), 
        runge_linear.(x_fine), 
        runge_spline.(x_fine)],
    linewidth=3)
plot(x_fine, [Runge.(x_fine), runge_linear.(x_fine)])
plot(x_fine, [Runge.(x_fine), runge_spline.(x_fine)])

################################################################################
# Interpolation 
################################################################################

####Linearly Interpolating log function
x = 0.1:2:10.1
y = log.(x)
x_fine = collect(0.1:0.1:10.1)
y_fine = log.(x_fine)

#plot
plot(x_fine, y_fine)
scatter!(x, y, markersize = 4)

#Linear Interpolation
interp_y = LinearInterpolation(x, y)
y_fine_interp = interp_y.(x_fine)
plot(x_fine, [y_fine, y_fine_interp])
scatter!(x, y, markersize = 4)

x2=[0.1,0.5,1.0,1.5,2.0,5.0,10.1]
y2 = log.(x2)
interp_y2 = LinearInterpolation(x2,y2)
y_interp = interp_y2.(x_fine)
plot(x_fine, [y_fine, interp_y2.(x_fine)])
##Extrapolation
log(11)
interp_y(11) #This will give you an error because we haven't allowed for extrapolation

#Need to set the extrapolation boundary coundition
#Line() will just extend the line from the last region
#Flat() will set everything equal to the last gridpoint. 
interp_y_extra = LinearInterpolation(x,y,extrapolation_bc = Line()) #Flat()
interp_y_extra(11)

x2=[0.1,0.5,1.0,1.5,2.0,5.0,10.1,19.0]
y2 = log.(x2)
interp_y2 = LinearInterpolation(x2,y2)
y_interp = interp_y2.(x_fine)
x_fine_2  = collect(0.1:0.1:20.1)
y_fine_2 = log.(x_fine_2)
plot(x_fine_2, y_fine_2)
scatter!([10.0,20.0], [log.(10.0),interp_y_extra(20)])
#Extrapolation works okay when close to actual gird, but can struggle far outside grid
log(11)
interp_y_extra(11)
log(20)
interp_y_extra(20)

#We can make our approximation better with uneven grid
x_uneven = 0.1 .+ 10.0 .*collect(0.0:0.2:1.0).^2
y_uneven = log.(x_uneven)
plot(x_fine, y_fine)
scatter!(x_uneven, y_uneven, markersize = 4)

interp_y = LinearInterpolation(x_uneven, y_uneven)
y_fine_interp = interp_y.(x_fine)
plot(x_fine, [y_fine, y_fine_interp])
scatter!(x_uneven, y_uneven, markersize = 4)


####Cubic Interpolating log function
x = 0.1:2:10.1
y = log.(x)
x_fine = collect(0.1:0.1:10.1)
y_fine = log.(x_fine)

#plot
plot(x_fine, y_fine)
scatter!(x, y, markersize = 4)

#Cubic Interpolation
interp_y = CubicSplineInterpolation(x, y)
y_fine_interp = interp_y.(x_fine)
plot(x_fine, [y_fine, y_fine_interp])
scatter!(x, y, markersize = 4)

###Extrapolation 
interp_y = CubicSplineInterpolation(x, y, extrapolation_bc = Line())
log(12)
interp_y(12)

#CubicSplineInterpolation does not support uneven grid in Interpolations.jl
x_uneven = 0.1 .+ 10.0 .*collect(0.0:0.2:1.0).^2
y_uneven = log.(x_uneven)
interp_y = CubicSplineInterpolation(x_uneven, y_uneven) #Error because x is not evenly spaced. 


################################################################
# Bilinear interpolation
################################################################

f(x, y) = 1 + x^2 + y^2 

grid_coarse = collect(0.0:1.0:5.0)
grid_fine = collect(0.0:0.01:5.0)
ncoarse, nfine = length(grid_coarse), length(grid_fine)
z_fine = zeros(nfine, nfine)

for i in eachindex(grid_fine), j in eachindex(grid_fine)
    x, y = grid_fine[i], grid_fine[j]
    z_fine[i,j] = f(x, y)
end

contourf(grid_fine, grid_fine, z_fine)

#now get to bilinear interp
z_coarse = zeros(ncoarse, ncoarse)
for i = 1:ncoarse, j = 1:ncoarse
    x, y = grid_coarse[i], grid_coarse[j]
    z_coarse[i,j] = f(x, y)
end

#create the interpolation
interp_z = LinearInterpolation((grid_coarse, grid_coarse), z_coarse) 
z_interp = zeros(nfine, nfine)
for i = eachindex(grid_fine), j = eachindex(grid_fine)
    z_interp[i,j] = interp_z(grid_fine[i], grid_fine[j])
end

#compare contour plots
contourf(grid_fine, grid_fine, z_interp) #slightly more jagged, but otherwise pretty good!
contourf(grid_fine, grid_fine, (z_interp - z_fine)) #Difference
contourf(grid_fine, grid_fine, (z_interp - z_fine)./z_fine) #Percent Difference
surface(grid_fine,grid_fine,z_fine)
surface!(grid_fine,grid_fine, z_interp)

################################################################################
#####Better version of optimal growth
######Optimal Growth with Interpolation / Optimization
################################################################################
###################################################################################################
# Optimal Investment Problem
###################################################################################################

using Parameters 
#This package gives us the @with_kw macro, which allows us to define default values in our structs

### Struct for our model paramters
@with_kw struct ModelParameters

    β::Float64 = 0.99 #We can define default values because of the @with_kw macro. 
    δ::Float64 = 0.025
    α::Float64 = 0.36

    k_grid::Vector{Float64}=collect(range(0.1, length = 100, stop = 45.0)) #### Lower the Number of grid points to 100      
    N_k::Int64 = length(k_grid)

    tol::Float64 = 10^-4

end


@with_kw struct ModelSolutions

    V::Vector{Float64}
    kp::Vector{Float64}

end 

function build_ModelSolutions(para)

    V = zeros(para.N_k)
    kp = zeros(para.N_k)

    sols = ModelSolutions(V,kp)

    return sols

end

function build_structs()

    para = ModelParameters()
    sols = build_ModelSolutions(para)

    return para, sols

end

### Bellman operator
function bellman(para, sols)
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    V_next = zeros(para.N_k)
    kp_next = zeros(N_k)

    #Interpolate value function for continuation value
    V_interp = LinearInterpolation(k_grid, V, extrapolation_bc = Line())

    for i_k = eachindex(k_grid)
        k = k_grid[i_k]
        budget = k^α + (1-δ)*k

        #Replace grid search with box constrained optimization!
        opt = optimize(kp -> -log(budget-kp) - β*V_interp(kp), 0.0, budget)

        V_next[i_k] = -opt.minimum
        kp_next[i_k] = opt.minimizer

    end

    return V_next, kp_next

end



### Solve model
function solve_model(para, sols)
    para, sols = build_structs();    
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    V_next = zeros(N_k)
    max_diff = tol + 10.0
    n = 0
    while max_diff > tol
        n +=1
        V_next, kp_next = bellman(para,sols)
        
        max_diff = maximum(abs.(V_next - V))
        sols.V .= V_next
        sols.kp .= kp_next

        @show n, max_diff

    end

    sols

end

para, sols = build_structs();

#So fast and just as accurate!
@time sols = solve_model(para,sols)
 

using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))