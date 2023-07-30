
########################################################################
# 0. Commments
########################################################################

# This is a commment. Julia will not run this line of code. 
# You can use comments to talk about what your code is doing. 
# Comment are an important part of having readable code. 

#=
This
is
a 
multiline
comment.
=#



########################################################################
# 1. Math
########################################################################

#Julia is just a fancy calculator.

2+2

10 + (20 - 24^3) / 5000 * 123



########################################################################
# 2. Suppressing and Printing Output
########################################################################

#This prints to the terminal
2+2

#Semi-colon will suppress output
2+3;
ans

#The print function will write output into REPL
print(2+2);

# The @show macro will print a line of code and its output to terminal. 
@show 2+2



########################################################################
# 3. Variables
########################################################################

# We can save the results of our calculations as variables
x = 2^3

y = x + 5

🐟 = 112

# There are special operations for updating variables
x += 1 #This adds on to x
x -= 1
x *= 2
x /= 10



########################################################################
# 4. Arrays, Vectors, and Matrices
########################################################################

# An array is a collection of elements. 

# A vector is a 1 dimensional array. 
v1 = rand(4) # Four random numbers between 0 and 1
v2 = collect(1:2:41)
v3 = collect(range(start = 1, stop = 10, length = 100))

# A matrix is 2-D array
A3 = rand(4,5)


# Julia has built-in matrix operations
x = rand(4)
Z = rand(4,4)

# Matrix multiplication
Z * x
x' * Z

#Transpose
Z'

#Matrix inverse
inv(Z)


# A . broadcasts a function to all the elements of an array
log.(Z) #Take log of all the elements of Z

########################################################################
# 5. Packages
########################################################################

# The right square bracket ] opens the package interface to add packages.
using Distributions, Plots



########################################################################
# 6. For Loops
########################################################################

# Construct data x_{t+1} = x_{t} + ϵ where ϵ ~ N(0,1)

x_vec = zeros(100)

for i in 2:100
    x_vec[i] = x_vec[i-1] + rand(Normal())
end

x_vec[4]

plot(x_vec)



########################################################################
# 6. Functions
########################################################################

# You can write your own functions

function OLS(Y,X)

     inv(X'*X)*(X'*Y)

end


coefs = zeros(100)

for i in 1:100
    ϵ_vec = rand(Normal(), 100)
    X = rand(Normal(), 100)
    Y = X + ϵ_vec

    coefs[i] = OLS(Y,X)
end

histogram(coefs)


########################################################################
# Booleans and Conditionals
########################################################################

# Booleans are true and false
true
false

2 == 3
5 > 1
6 <= 7
10 >= 100
7 != 6


# Conditional statements

x = 3

if x > 5
    print("Yes")
else
    print("No")
end


x = 0
while x < 100
    x += 1
end
@show x













































