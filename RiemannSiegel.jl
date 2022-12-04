#=**
**
**    Riemann-Siegel formula for roots of Zeta(s) on the critical line.
**
**************************************************************************
**    Matthew Kehoe
**    10/30/2022
**
**    This program finds roots of Zeta(s) using the Riemann-Siegel
**    formula. The Riemann–Siegel theta function is approximated by
**    Stirling's approximation. It also uses an interpolation method to
**    locate zeros. The coefficients for R(t) are handled by the Taylor
**    Series approximation originally listed by Haselgrove in the 1960s.
**    It is necessary to approximate these coefficients in order to
**    increase computational speed.

Future Work

1. Implement a cleaner version of findRoots(). Investigate interval methods.
2. Use GPU and parallel processing capability. A lot of good ideas are discussed in
https://www.youtube.com/watch?v=8E0qdO_jRZ0.
3. Calculate all of the Riemann-Siegel coefficients explicitly.
4. Implement the Odlyzko–Schönhage algorithm.
**************************************************************************=#

using Roots
using Printf
# using Quadmath - for Float128. Not currently implemented.

# Use the functions defined in zetaFunctions.jl
include("zetaFunctions.jl")

"""
Search for zeros of the Riemann-Siegel Z function in the interval [a,b] by
the find_zeros function in Roots.jl.
### Input
    @param a - the starting value of the interval.
    @param b - the ending value of the interval.
### Output
    @return an array of zeros of zeta(s) on the critical line.
"""
function findRoots(a::Float64, b::Float64, intervalSize::Float64)
    # This function is inaccurate when the interval size is large
    startPoint = a
    endPoint = (b % intervalSize) + intervalSize
    numIntervals = b ÷ intervalSize
    rtsArray = zeros(0)
    @inbounds for i=1:numIntervals
        f(x) = RiemennZ(x, 4);
        tempArray = find_zeros(f, startPoint, endPoint)
        startPoint = endPoint
        endPoint = endPoint + intervalSize
        append!(rtsArray, tempArray)
    end
    for root in rtsArray
        println(root)
    end
    println(size(rtsArray))
end


"""
    Extremely inefficient way of calculating roots manually.
"""
function findRootsSlow(a::Float64, b::Float64)
    f(x) = RiemennZ(x, 4);
    count = 0
    step_size = 0.01
    @inbounds while a < b
        a += step_size
        # This is not the most efficient way of looking for roots
        if sign(f(a)) != sign(f(a + step_size))
            x = find_zero(f, [a, a + step_size], Order5())
            # Print the roots to the console
            println(x)
            # Save roots to a file
            #io = open("roots.txt","a")
            #    println(io,x)
            #close(io)
            count = count + 1
        end
    end
    @printf("Number of roots in interval: [%d,%d] is %d\n", 0, b, count)
end


# Main
function main()
    println("Zeroes on the critical line Zeta(1/2 + it).")
    a = 1.0
    b = 100.0
    #b = 74920.9                 # Should return 100,000 zeros!!!
    intervalSize = 1.0        # Setting intervalSize = 1.0 returns 100,000 zeros correctly
    @time findRootsSlow(a, b)
    #@time findRoots(a, b, intervalSize)
    
end

main()
