#=/**************************************************************************
**
**    Bernoulli numbers to calculate Zeta(2n) and Zeta(-n)
**
**************************************************************************
**    Matthew Kehoe
**    11/26/2022
**
**    The Bernoulli numbers are related to functional values of the Riemann
**    zeta function. One of these relationships, originally found by Euler,
**    states that Zeta(2n) = -1^(n+1) * B_2n * (2*pi)^2n / (2 * (2n)!) for
**    positive integers n. Another relationship for negative integer values
**    of n is that Zeta(-n) = -B_n+1 / (n+1). This program handles both cases
**    by calculating the nth Bernoulli number through the Akiyama-Tanigawa
**    algorithm and printing out the results to the screen.
**************************************************************************/=#

"""
Rational representation of the Zeta function for positive even integers. 
### Input
    * @param n the integer value of Zeta(2n)
    * Numerator = (-1)^(n+1) * B_2n * (2pi)^(2n)
    * Denominator = 2*(2n)!
### Output
    * @return the numerator of Zeta(2n)
    * @return the denominator of Zeta(2n)
    * @return the value of Zeta(2n)
"""
function evenZeta(n) 
    top = (BigInt(-1))^(n+1) * Bernoulli(2*n) * (BigInt(2))^(2*n)
    bottom = 2 * factorial(big(2 * n))
    val = top / bottom
    num = numerator(val)
    den = denominator(val)
    return num, den, val
end


"""
Rational representation of the Zeta function for negative integers. 
### Inputs
    * @param n the integer value of Zeta(-n)
    * Numerator = -B_n+1
    * Denominator = n+1
### Outputs
    * @return the numerator of Zeta(-n)
    * @return the denominator of Zeta(-n)
    * @return the value of Zeta(-n)
"""
function oddZeta(n)
    val = -Bernoulli(n+1) / (n+1)
    num = numerator(val)
    den = denominator(val)
    return num, den, val
end

"""
Compute the  Bernoulli number of an integer through the Akiyama-Tanigawa algorithm. 
### Inputs
    * @param n the integer
### Outputs
* @return the Bernoulli number, B_n
"""
function Bernoulli(n)
    A = Vector{Rational{BigInt}}(undef, n + 1)
    for m = 0 : n
        # use double divides to store the result as a rational integer
        A[m + 1] = 1 // (m + 1)
        for j = m : -1 : 1
            A[j] = j * (A[j] - A[j + 1])
        end
    end
    return A[1]
end

""" 
Print out the values of Zeta(2n) for even positive integers.
### Inputs
    * @param n the integer
### Outputs
* @return all positive even integer values of Zeta(2n) in rational form
"""
function printZeta2N(n)
    for i in range(0, n)
        num,den,val = evenZeta(i)
        if abs(num) != 1
            println("ζ(" * string(2*i) * ") = " * string(num) * string(π) * 
                    "^" * string(2*i) * " / " * string(den))
        elseif i == 0
            println("ζ(" * string(2*i) * ") = -" * string(π) * 
                    "^" * string(2*i) * " / " * string(den))
        else
            println("ζ(" * string(2*i) * ") = " * string(π) * 
                    "^" * string(2*i) * " / " * string(den))
        end
    end
end

""" 
Print out the values of Zeta(-n) for odd integers.
### Inputs
    * @param n the integer
### Outputs
* @return all negative values of Zeta(-n) in rational form
"""
function printZetaOddN(n)
    for i in range(0, n)
        num,den,val = oddZeta(i)
        if num % 2 == 0
            println("ζ(" * string(-i) * ") = " * string(0))
        else 
            println("ζ(" * string(-i) * ") = " * string(num) 
                            * " / " * string(den))
        end
    end
end

# Main - Print out values of Zeta(k) up to the desired integer value of k.
function main()
    n = 20
    println("Even integer values of ζ(2n) up to n = " * string(n) * ".")
    println()
    printZeta2N(n)
    println()
    println("Negative integer values of ζ(-n) up to n = " * string(n) * ".")
    println()
    printZetaOddN(n)
end

main()