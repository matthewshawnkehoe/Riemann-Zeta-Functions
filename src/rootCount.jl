#=**
**
**    Count the number of roots of Z(t) by the number of Gram points.
**
**************************************************************************
**    Matthew Kehoe
**    12/04/2022
**
**    TODO: Rewrite the paragraph below
**    This program counts the number of zeros of the Riemann zeta function 
**    between Graham points BEGIN and END. The program establishes the         
**    existence of a zero with complete certainty, and thus produces a lower   
**    bound on the number of zeros in the search range. There is the           
**    possibility that a zero will be missed in the case where |Z(t)| is       
**    less than the error tolerance, and hence the sign of Z(t) cannot be     
**    determined with certainty. In these situations, a warning is written     
**    to the RESULTFILE indicating that the questionable region should be      
**    examined with more accuracy.  
**************************************************************************/=#


"""
Locate the n^th Gram point through Newton's method.
### Input
@param n - the n^th Gram point.
### Output
@return the value of the n^th Gram point.
"""
function gram(n::Int)
  tn = 0.0
  tn1 = 0.5*n + 20.0

  while abs(tn - tn1) > 1e-12
      tn = tn1
      # Compute t_1 = t_0 - f(t_0) / f'(t_0) 
      # Here, f is the Riemann-Siegel theta function minus n times pi
      tn1 = tn - (theta(tn) - n*pi) / thetaDeriv(tn)
  end

  return tn1
end


"""
    Riemann-Siegel theta function using the approximation by the
    Stirling series.
### Input
    @param t - the value of t inside the Z(t) function.
### Output
    @return Stirling's approximation for theta(t).
"""
function theta(t::Float64)
    twopi = 2.0 * pi
    return (t/2.0 * log(t/twopi) - t/2.0 - pi/8.0
            + 1.0/(48.0*t) + 7.0/(5760*t^3))
end


"""
    Derivative of Riemann-Siegel theta function using the approximation by the
    Stirling series.
### Input
    @param t - the value of t inside the Z(t) function.
### Output
    @return Derivative of Stirling's approximation for theta(t).
"""
function thetaDeriv(t::Float64)
    twopi = 2.0 * pi
    return (0.5 * log(t/twopi) - 1.0/(48.0 * t^2) - 7.0/(1920.0*t^4))
end


"""
Find the expected number of roots between Graham points n and m. If we miss a root
in a Gram block, then we will write out error information to an output file.
### Input
@param n - the n^th Gram point.
@param m - the m^th Gram point.
@param outFile - an error file describing where we might have missed a root.
### Output
@return the number of roots between Gram points n and m.
"""
function gramBlock(n::Int64, m::Int64, outFile)

    #int num;                        /* root counter                         */
    #int p;                          /* partitions between graham points     */
    #int j;                          /* index of current graham point        */
    #int k;                          /* index of current partition point     */

    #long double Z(long double,int); /* Riemann-Siegel Z function            */
    #long double graham(int);        /* Graham point evaluation function     */
    #long double gp1;                /* Graham interval partition point      */
    #long double gp2;                /* Graham interval partition point      */
    p = 1
    num = 0

    # If the number of roots is less than the expected number of roots
    while (num < (m - n) && p < 1024) 

        # Redfine the Gram point partition and start at the beginning of the 
        # Gram block
        num = 0                                    
        p = 2*p                                   
        j = n;                                      
        gp2 = gram(j)   
        
        # For every Gram point besides the last
        while (j <= (m-1)) 

          # Step through the successive partition points
          k = 0
          while (k <= (p-1))  
            
            # Determine coordinates of the adjacent partition points
            gp1 = gp2

            gp2 = gram(j) * (1 - (k+1)) / p + gram(j+1) * (1 - (k+1)) / p

            # Evaluate Z(t) at adjacent partiton points
            
            
          end

        end
  

    end

end

function main()
    println("Gram points")
    for i in -1:3
        println(gram(i))
    end
end

main()