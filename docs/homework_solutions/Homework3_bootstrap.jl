
################################################################################
# Bootstrap
################################################################################

using Distributed
addprocs(10) # Add 10 workers

@everywhere using CSV
@everywhere using DataFrames 
@everywhere using GLM #Regressions
@everywhere using Optim
@everywhere using Random
@everywhere using SharedArrays
@everywhere using StatsBase

#Load Data set on all cores
@everywhere df = CSV.read("lwage.csv", DataFrame)
@everywhere df.exp_sq = df.exp .^2

#Define Log Likelihood on all cores
@everywhere function loglikelihood_LR(n, X, Y, θ)

    #n = length(Y)

    #Unpack θ
    #β_0 = θ[1]
    #β_1 = θ[2]
    #β_2 = θ[3]
    #β_3 = θ[4]
    #σ = θ[5]

    #Log Likelihood
    -n/2*log(2π) - n/2* log(θ[5]^2) - (sum((Y .- θ[1] - X * θ[2:4] ).^2) / (2*θ[5]^2)) 

end

# Run main estimation on full data 
lmsum = lm(@formula(lwage ~ educ + exp + exp_sq), df)
coef_OLS = coef(lmsum)
σ_OLS = std(residuals(lmsum))
Y = Vector{Float64}(df.lwage)
X = Matrix{Float64}(df[:,2:4])
opt = optimize(θ -> -loglikelihood_LR(length(Y),X, Y, θ), [coef_OLS; σ_OLS])

θ_hat = opt.minimizer


# Define bootstrap procedure on all cores
@everywhere function boot_once(df, init_guess)

    n = nrow(df)
    df_b = df[sample(1:n,n),:]

    Y_b = Vector{Float64}(df_b.lwage)
    X_b = Matrix{Float64}(df_b[:,2:4])

    opt = optimize(θ -> -loglikelihood_LR(length(Y_b),X_b, Y_b, θ), init_guess)

    opt.minimizer

end


function boot(N_b, θ_hat)

    results = SharedArray{Float64}(zeros(N_b, length(θ_hat)))

    @sync @distributed for i = 1:N_b
        Random.seed!(20230813 + i)
        results[i,:] = boot_once(df,θ_hat)
    end

    results

end

results = boot(100,θ_hat)


@show ste_boot = [std(results[:,i]) for i in 1:5]








