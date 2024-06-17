
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

#emojis!!
🐟 = 112
🤑 = 300
🐟 + 🤑
#type '\', then enter the desired emoji name

#you can also do greek letters:
β= 3

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

#if you want to create your own (column) vector with specified values, just enter it like this:
v1 = [4.0, 0.3, 17.0, -0.2]
#note the comma between the entries
#the comma tells julia that the next entry should be in the next row 

#same, but for a row vector
v1 = [4.0 0.3 17.0 -0.2]
#note the lack of commas or any separator

# A matrix is 2-D array
A3 = rand(4,5)

#to create one of these manually, you can do it like this:
A4 = [1.0 2.0; 3.0 4.0]
#(this is my preferred way)

#OR you could do it like this: (needlessly confusing, in my opinion)
A4 = [[1.0, 3.0]  [2.0, 4.0]]

#OR you could even do this (please don't do this, even though it technically works it just feels wrong)
A5 = [1.0 2.0
3.0 4.0]

#the confusing thing for me to learn initially: knowing which way I should enter rows vs. columns - I would ~highly~ recommend printing the matrix you just created to the console to make sure that the matrix actually looks like what you want it to look like

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

# The right square bracket ] opens the package interface to add packages. To install a package, you will need to type the ] (which will then show the prompt 'pkg>' instead of 'julia>'. From there, you will want to type add 'package_name', where 'package_name' is the name of your package. Once installed, you can load the package by typing 'using package_name')

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

#by the way: the first time you plot something, it takes longer than subsequent times
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
    ε_vec = rand(Normal(), 100)
    X = rand(Normal(), 100)
    Y = X + ε_vec

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













































