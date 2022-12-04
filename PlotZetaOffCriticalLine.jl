#=/**************************************************************************
**
**    Plot real values of the Zeta Function  
**
**************************************************************************
**    Matthew Kehoe
**    11/12/2022
**
**    This program plots the value of Zeta(s) for real s. Different values
**    of s are computed through truncating the infinite series representation, 
**    Riemann's functional equation, and the Lanczos and Stirling approximation
**    of the Gamma function. he current implementation does not support large
**    negative real values of s.
**************************************************************************/=#

using Plots
using LaTeXStrings
using Printf

"""
The standard truncation of an infinite series for the Zeta function. The
summation stops after 10,000 iterations.
### Input
	 @param s - the value of s.
### Output
	 @return the value of Zeta(s) through an approximation of the infinite series.
"""
function standardZeta(s, limit=10000)
    result = 0
    for n in 1:limit
        result += (1/n)^s
    end
    return result
end


"""
The Dirichlet series approximation of the Zeta function which is valid for Re(s) > -1. 
The maximum number of iterations is set to 10,000.
### Input
	 @param s - the value of s.
### Output
	 @return the value of Zeta(s) through the Dirichlet series approximation
"""
function modifiedZeta(s, limit=10000)
    cons = 1 / (s - 1)
    result = 0
    for n in 1:limit
        firstTerm = (n * (n + 1)) / 2
        secondTerm = (2*n + 3 + s) / (n + 1)^(s+2)
        thirdTerm = (2*n-1-s) / n^(s+2)
        result += firstTerm * (secondTerm - thirdTerm)
    end
    return cons * result
end

"""
Approximation of the Gamma function by the Lanczos Approximation.
### Input
    @param s - the value of s.
### Output
    @return the approximate value of Gamma(s).
"""

function lancGamma(s)

    # Coefficents as calculated by Godfrey, Paul (2001). "Lanczos implementation of the 
    # gamma function". Numericana
    p = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012,
        9.9843695780195716e-6, 1.5056327351493116e-7]

    s = complex(s)

    # Remove the imaginary part if it is small
    EPSILON = 1e-07
    function drop_imag(s)
        if abs(imag(s)) <= EPSILON
            s = real(s)
        end
        return s
    end

    # If real(s) < 0.5, use the reflection formula. Otherwise, approximate the gamma 
    # function through the lanczos approximation.
    if real(s) < 0.5
        gamma = pi / (sin(pi * s) * gamma(1 - s))  # Euler's Reflection formula
    else
        s -= 1
        Ag = 0.99999999999980993
        for (i, pval) in enumerate(p)
            Ag += pval / (s + i)
        end
        t = s + length(p) - 0.5
        gamma = sqrt(2 * pi) * (t ^ (s + 0.5)) * exp(-t) * Ag
    end
    return drop_imag(gamma)
end

"""
Approximation of the Gamma function by the Stirling Approximation. 
### Input
    @param s - the value of s.
### Output
	@return the approximate value of Gamma(s).
"""
function strlGamma(s::Float64)
    s = s+1
    stirlingApprox = sqrt(2 * pi/s) * (s/exp(1) * sqrt(s * sinh(1/s) + 1/(810*s^6))) ^ s
    return stirlingApprox
end


"""
Calculate Zeta(s) for real s.

1. The first if statement handles when s < 0 and s is a multiple
 of 2k. These are trivial zeroes where Zeta(s) = 0.

2. The second and third if statements handle the case where s = -1 
or s = 1. In the later case, Zeta(s) is undefined. 

3. The fourth if statement handles when s < -170. The value of Zeta(s)
is calculated by recursion through the functional equation and
Stirling's approximation for Gamma(s).

4. The fifth if statement handles when s < -1. The value of Zeta(s)
 is calculated by recursion through the functional equation and
 the Lanczos approximation for Gamma(s).

5. The last if statement calculates Zeta(s) through the 
 standard truncation of the infinite sum provided s > 1.
 
### Input
* @param s - the value of s.
### Output
* @return the value of Zeta(s)
"""
function ζ(s)
    term = 2^s * pi^(s-1) * sin((pi*s)/2)

    if (s <= -2 && s % 2 == 0)
        return 0
    elseif s == 1
        return Inf
    elseif s == -1
        return -1/12
    elseif (s < -170)
        return term*standardZeta(1-s)*strlGamma(1-s)
    elseif (s < -1)
        return term*standardZeta(1-s)*lancGamma(1-s)
    elseif (s > -1 && s < 1)
        return modifiedZeta(s)
    else 
        return standardZeta(s)
    end
end
    
function main()
    a = -30         # lower bound of real domain
    b = -28          # upper bound of real domain
    numSamples = 100     # number of samples to compute

    samples = numSamples - 1
    sampleIndex = 0 
    intervalLength = b - a

    rePart = Float64[]
    axisVal = Float64[]
    while sampleIndex <= samples
        s = a + sampleIndex/samples * intervalLength
        zeta = ζ(s)

        #@printf("s = %.12f\t zeta(s) = %.12f \n",s,zeta)

        append!(rePart, zeta)
        append!(axisVal, s)
        
        sampleIndex = sampleIndex + 1
    end

    plot((axisVal,rePart), color=:blue, width=3, title=L"Re(\zeta(s))", xaxis=L"s", 
          yaxis = L"\zeta(s)", framestyle= :origin,label=L"\zeta(s)")
    #scatter!((rePart),color=:red, label="Points")
end

main()

