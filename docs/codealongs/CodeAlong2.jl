
###################################################################################################
# Machine code
###################################################################################################

# Julia is the translator between you and your computer. You learn the Julia programming language, 
# and Julia converts it into machine language that your computer understands. If you ever want to 
# see what Julia is say to your computer, you can use the @code_native macro. 

@code_native 2+2

@code_native (1 + 8*2)^2

# That garbled nonsense that appears in the terminal is machine language.



###################################################################################################
# Local Scope vs. Global Scope
###################################################################################################

# The scope of a variable is where in the code the region is accessible. 

# Global scope are essentially anything that are defined outside of 
# functions, for loops, conditionals, etc. 

x = 2 #Here, x is a global variable. 

# We can operate on x in the global scope. 
x + 4

# When you make a function, you define a new local scope.
function f()
    # Things defined inside this function are only accessible within our function. 
    z = 3
    return z
end

f()

# Inside the local scope of the function f, we defined the variable 'z'. 
# But we can't access 'z' outside of the function in the global scope
z 
# You get an error if you try to use z in global scope because it only exists in the local scope of 
# the function f. 


#For loops also define a local scope. 
for i in 1:10
    y = 7
    @show y + i 
end

#y is not defined in global scope. Only in the local scope of the above for loop.
y

# If we define x in global scope, then it will be updated when we operated on it in local scope. 
x

for i in 1:10
    x = 7
    @show x + i 
end

# x is now 7
x




###################################################################################################
# Tip 1: Write things inside functions
###################################################################################################

# Julia is fast because it performs just-in-time compilation. This means when you run your code, 
# Julia will compile it and then run it. A piece of code only needs to be compiled once. This is why
# things in Julia will often be slow the first time you run them, but much faster the second time. 
# The first run includes both the compilation and the execution. The second time, you are only
# executing the code. 



###Example 1: Very bad: Everything here is defined in global. 
# We are not taking advantage of compilation at all here. 

m = 3
count = 0

@time for i = 1:10^7
    count += m*i
    count -= m*i
end


###Example 2: A little better

function update(x,m,i)
    x += m*i
    x -= m*i
    return x
end

@time for i = 1:10^7
    count = update(count,m,i)
end



###Example 3: Good, used a function. Everything happens in local scope
function g(m)

    count = 0.0

    for i = 1:10^7
        count = update(count, m, i)
    end

    return count

end


@time g(m)



###################################################################################################
# Types and Structs
###################################################################################################

# Types describes the structure and behavior of an element. 
# We have seen lots of types already in the first  lecture

# Integer with 64 bits
typeof(1)
# Floating number with 64 bits
typeof(1.0)
# String 
typeof("abc")
# Character
typeof('a')
# Vector of Floats, which is a 1-D array
typeof(zeros(2))
# Matrix of Floats, which is a 2-D Array
typeof(ones(3,2))
# 4-D Array of integers
typeof(zeros(Int64, 3,3,3,3))


#You can define the type of a variable using "::" when you call it. 
x = 32.0::Float64 
typeof(x)
#This is entirely unneccesary, though, because Julia will figure out your types for you. 
y = 32.0
typeof(y)

# You can also specify the type of your function arguments.
function h(x::Int64)

    return x * 100

end

h(10) #Your function works with integers
h(10.0) #But it will not work with floats. 

# You can avoid this by allowing julia to infer the type

function h2(x)
    return x * 100
end

#Now this will work with both types. 
h2(10)
h2(10.0)

# You can define your own types using the struct ("structure") keyword
# This type has two fields, a and b. 
struct MyType
    a 
    b
end

#We can define an instance of MyType
example = MyType(25.0, "abc")

#We can access the fields of example like this
example.a
example.b

# Structs are very very useful when coding up big models with many parameters or output. 
# Instead of carrying around all these different things, you just need to carry around your struct 
# container.

###################################################################################################
# Tip 2: Declare concrete types inside your structs. 
###################################################################################################

# Let's define a struct here without declaring the type of the element. 
# Here, "a" could be anthying. A float, integer, string, character, matrix, or anything at all.  
struct Container_Any
    a
end

# Here, we specified that "a" must be a float. If we don't make assign a float to "a", we will get
# an error. 
struct Container_Float
    a::Float64
end


#Let's assign a float to both anyways. 
c_any = Container_Any(1.0)
c_float = Container_Float(1.0)
typeof(c_any)
typeof(c_float)
typeof(c_any.a)
typeof(c_float.a)

#Let's define a function that operates on our containers. 
function f(x)
    tmp = x.a
    #A bunch of calculations
    for i in 1:10^7
        tmp = i + x.a
    end
    return tmp
end

#Let's see how fast our function is
@time f(c_any)
@time f(c_float)

#But now run them a second time
@time f(c_any)
@time f(c_float)
# Now f is much faster on c2 than on c1
# This is because f has now been compiled. Because the type is declared for element within c2, 
# Julia is able to write very efficient machine code to execute the program instantly. 
# While for c1, it doesn't know the type will always be float, and has to do many useless calculations. 



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

    k_grid::Vector{Float64}=collect(range(0.1, length = 1800, stop = 45.0))       
    N_k::Int64 = length(k_grid)

    tol::Float64 = 10^-4

end


struct ModelSolutions

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

    for i_k = eachindex(k_grid)
        max_util = -1e10
        k = k_grid[i_k]
        budget = k^α + (1-δ)*k

        for i_kp = eachindex(k_grid)
            
            c = budget - k_grid[i_kp]

            if c > 0
                
                V_temp = log(c) + β*V[i_kp]
            
                if V_temp > max_util
                    max_util = V_temp
                    kp_next[i_k] = k_grid[i_kp]
                end
                V_next[i_k] = max_util
            end

        end

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
        V .= V_next
        kp .= kp_next

        @show n, max_diff

    end

end

para, sols = build_structs();

@time solve_model(para,sols)
 

using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))


