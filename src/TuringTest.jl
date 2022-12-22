#=**
**
**    Apply Turing's method to bound S(g_n) by calculating the adjustments in
**    h_n. 
**
**************************************************************************
**    Matthew Kehoe
**    12/22/2022
**
**    Program to find the adjustments h_m used in Turing's method for bounding 
**    S(g_n)= N(g_n) - theta(g_n) / pi - 1.
**
**    For the final Gram point, g_END, the program finds values of h_(END-SAMPLES+1),      
**    ..., h_END such that 
**
**    1.)  g_(END-SAMPLES+1)+h_(END-SAMPLES+1), ..., g_END+h_END is            
**         strictly increasing, and 
**
**    2.) (-1)^n * Z(g_n + h_n) > 0 for n = END-SAMPLES+1, ..., END
**
**    Output the results to a file consisting of n, g_n, h_n, g_n + h_n, Z(g_n + h_n),   
**    for n = END-SAMPLES+1, ..., END. 
**************************************************************************/=#

# Use the functions defined in zetaFunctions.jl
include("zetaFunctions.jl")

using Printf


function main()
    g_END = 12193892             # Final Gram point 
    SAMPLES = 20                 # Number of samples to compute

    g_n = Array{Float64,1}(UndefInitializer(), SAMPLES)   # Array of Gram points
    h_n = Array{Float64,1}(UndefInitializer(), SAMPLES)   # Array of adjustment terms

    # Initialize Variables
    num_samples = SAMPLES                          
    j = 0                                            
    k = 0                                         
    m = 0                                         
    g_initial = g_END - num_samples + 1             
    step      = 0.1                          
    h_total   = 0.0
    
    info_file = open("turing.txt", "w")

    while (j < num_samples)
        # Evaluate Gram points in range and initialize adjustment vector
        g_n[j+1] = gram(g_initial + j)
        h_n[j+1] = 0.0
        j += 1
    end

    while (k < num_samples)
        # For every Gram point, set the number of steps to zero
        n = 0

        # The error tolerance is epsilon(t) <= 0.053t^(-5/4) for Gram points
        leftArg = g_n[k+1] - n*step
        rightArg = g_n[k+1] + n*step
        epsLeft = 0.053 / ((leftArg) ^ 1.25)
        epsRight = 0.053 / ((rightArg) ^ 1.25)

        # Keep on stepping one step to the left/right of g_n until (-1)^n * Z(g_n + h_n) > 0 
        # is satisfied where h_n := +/- n * step
        sign = even(g_initial+k)
        while (sign * RiemennZ(g_n[k+1]-n*step,4) < epsLeft  &&  
                        sign * RiemennZ(g_n[k+1]+n*step,4) < epsRight)
            n += 1
        end
        
        # Track the minimal adjustment found in the array of adjustment terms
        if (sign * RiemennZ(g_n[k+1]-n*step,4) > epsLeft)
            h_n[k+1] = -n * step
        else
            h_n[k+1] = n * step
        end

        # Add the individual adjustments to the h_total
        h_total += h_n[k+1]

        # If the sequence g_k + h_k isn't increasing, reduce the step size and start over.
        # Otherwise, go to the next Gram point.
        if ( k > 1 && (g_n[k]+h_n[k] <= g_n[k-1]+h_n[k-1]))
            h_total = 0
            step = step / 2
            k = 1
        else
            k += 1
        end
    end

    # Write out the results
    while (m < num_samples)
        @printf( info_file, "%d & %11.3Lf & %5.3Lf &  %11.3Lf & %5.3Lf \n", 
                 g_initial + m,
                 g_n[m+1],
                 h_n[m+1],
                 g_n[m+1]+h_n[m+1],
                 RiemennZ(g_n[m+1]+h_n[m+1],4) )
        m += 1
    end

    @printf(info_file, "h_total = %5.3Lf\n", h_total)

    # Close the output file
    close(info_file)

    @printf("%20.10Lf, %20.10Lf \n", gram(1181229), gram(1181235))
end

main()