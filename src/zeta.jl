# Simple program to calculate the value of zeta(s) when Re(s) > 1

using Plots

function ζ(z,limit=1000)
    result = 0
    for i in 1:limit
        result += (1/i)^z
    end
    return result
end

println(ζ(2,10_0000_000)) # Returns 1.644934057834575
# Basel Problem
println(pi^2/6) # Returns 1.6449340668482264

# Imaginary value
println(ζ(2+5im)) # Returns 0.8510045028264933 + 0.09880538410302253im

# Add complex numbers
function ζComplex(z,limit=1000;point_ar = false)
    points = ComplexF64[0]
    result = 0
    for i in 1:limit
        result += (1/i)^z
        if point_ar
            push!(points,result)
        end
    end
    return result, points
end

z = 0.5+14.1347251417346937904572519835625im

println(ζComplex(z,100000)) # Returns 1.644934057834575


test = 9