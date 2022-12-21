#=**
**
**    Verify the Riemann hypothesis up to a final Gram point. The Riemann-
**    Siegel formula is used for fast evaluation of the Z(t) function. 
**
**************************************************************************
**    Matthew Kehoe
**    12/04/2022
**
**    This program counts the number of zeros of the Riemann zeta function 
**    between Gram points START and END. In doing so, it establishes the        
**    existence of all zeros up to a specified end point and therefore produces 
**    a lower bound on the number of zeros in the search range. In the 
**    possibility that a zero is missed in the case where |Z(t)| is less than 
**    the error tolerance epsilon, a warning is written to an information file 
**    indicating that this region should be analyzed with more precision. The 
**    epsilon used in this program is that of eps(t) <= 0.053*t^(-5/4), as 
**    suggested by A. Odlyzko.
**************************************************************************/=#

# Use the functions defined in zetaFunctions.jl
include("zetaFunctions.jl")

using Printf


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

    # Variables

    # num: The number of roots between Gram points n and m.
    # p: The number of partitions between Gram points.
    # j: Index of the current Gram point.
    # k: Index of the current partition point.
    # Z(t,r): The Riemann-Siegel Z(t) function with r remainder terms
    # gram(n): Evaluation of the n^th Gram point.
    # gp1, gp2: Gram interval partition points
    
    num = 0 
    p = 1

    # If the number of roots is less than the expected number of roots
    while (num < (m - n) && p < 1024) 

        # Redefine the Gram point partition and start at the beginning of the Gram block
        num = 0                                    
        p = 2*p                                   
        j = n                                      
        gp2 = gram(j)   

        # Find the number of points in the Gram block [g_n, g_m) using the C0 and C1
        # remainder terms
        stage = 1
        num = numPointsInGramBlock(j, n, m, p, num, gp2, stage)
    end

    # If the number of roots is still less than the expected number of roots, reset the partition size
    p = 1

    # And start over using all the remainder terms for a more accurate Z(t)
    while (num < (m - n)) 

        # Once again redefine the Gram point partition and start at the beginning of the Gram block
        num = 0                                    
        p = 2*p                                   
        j = n                                      
        gp2 = gram(j) 
        
        # Find the remaining number of points in the the Gram block [g_n, g_m) using all remainder terms
        stage = 2
        num = numPointsInGramBlock(j, n, m, p, num, gp2, stage) 

        # Break if partition depth is too deep
        if (p == 65536)
            break
        end
    end                    

    # If the number of roots isn't the expected number, write a warning message to the output
    # file. Rosser's rule states that Gram blocks of length k should contain at least k zeros
    if (num < (m - n))
        @printf(outFile, "Partition depth reached in Gram block [%d,%d], possible 
                violation of Rosser's rule.\n",n,m)
    end
    
    return num
end


"""
Intermediate calculations to find the number of roots in the Gram block [g_n, g_m).
### Input
    @param n - the n^th Gram point.
    @param m - the m^th Gram point.
    @param p - the number of partitions between Gram points
    @param outFile - an error file describing where we might have missed a root.
    @param num - Current number of roots between Gram points n and m.
    @param stage - Decision to use two remainder terms or all remainder terms.
### Output
    @return the number of roots between Gram points n and m.
"""
function numPointsInGramBlock(j::Int64, n::Int64, m::Int64, p::Int64, num::Int64, gp2::Float64, stage::Int64)
        
    # For every Gram point besides the last
    while (j <= (m-1)) 

        # Step through the successive partition points
        k = 0
        while (k <= (p-1))  
            
            # Determine coordinates of the adjacent partition points
            gp1 = gp2
            gp2 = gram(j) * (1 - (k+1)) / p + gram(j+1) * (k+1) / p

            # The error tolerance is epsilon(t) <= 0.053t^(-5/4) for both Gram points
            eps1 = 0.053 / (gp1 ^ 1.25)
            eps2 = 0.053 / (gp2 ^ 1.25)

            # Break up the evaluation into two stages using only two remainder terms in 
            # the first stage. Then, evaluate Z(t) at adjacent partiton points and, if the signs differ,
            # increment the root count
            if stage == 1
                if (RiemennZ(gp1,1) > eps1 && RiemennZ(gp2,1) < -eps2)
                    num += 1
                elseif (RiemennZ(gp1,1) < -eps1 && RiemennZ(gp2,1) > eps2)
                    num += 1
                end
            else 
                if (RiemennZ(gp1,4) > eps1 && RiemennZ(gp2,4) < -eps2)
                    num += 1
                elseif (RiemennZ(gp1,4) < -eps1 && RiemennZ(gp2,4) > eps2)
                    num += 1
                end
            end

            # Then, go to the next partition point
            k += 1          
        end

        # After completing the above procedure, go to the next Gram point
        j += 1 
    end    

    return num
end


function main()

    START = -1    # Initial Gram point
    END   = 3000    # Final Gram point

    # Initialize Variables
    current_gpt            = START     # Index of the current Gram point                  
    final_gpt              = END       # Index of the final Gram point                 
    tot_num_roots          = 0         # Total number of roots found           
    gblock_roots           = 0         # Number of roots in the current Gram block                 
    tot_num_gblocks        = 0         # Total number of Gram blocks found                
    tot_num_gblock_roots   = 0         # Total number of roots found in the Gram blocks                
    largest_gblock_start   = -1        # Largest Gram block start                 
    largest_gblock_end     = -1        # Largest Gram block end
    
    info_file = open("gramInfo.txt", "w")

    # If we are not at the last Gram point
    while (current_gpt < final_gpt) 
        # Go to the next Gram point
        current_gpt += 1

        # The error tolerance is epsilon(t) <= 0.053t^(-5/4)
        eps = 0.053 / (gram(current_gpt) ^ 1.25)

        # If Grams Law is satisfied, increment the total number of roots
        if (even(current_gpt)*RiemennZ(gram(current_gpt),1) >  eps)
            tot_num_roots += 1
        # Otherwise, increase the accuracy of Z(t) by using more remainder terms
        elseif (even(current_gpt)*RiemennZ(gram(current_gpt),4) >  eps)
            tot_num_roots += 1
        # If both cases fail, we need to compute a Gram block.
        else
            # Find the left index of the Gram block
            gblock_start = current_gpt - 1

            # And test successive Gram points until Gram's law is satisfied
            while (even(current_gpt)*RiemennZ(gram(current_gpt),1) <=  eps)
                current_gpt += 1
            end

            # After this is complete, we should be at the right index of the Gram block
            gblock_end = current_gpt

            # And should be able to compute the bumber of roots in the Gram block
            gblock_roots = gramBlock(gblock_start, gblock_end, info_file)

            # Increment the root and Gram block counters
            tot_num_gblocks += 1
            tot_num_gblock_roots += gblock_roots
            tot_num_roots += gblock_roots

            # Also, determine which is the Gram block is the largest found so far
            if (gblock_end - gblock_start > largest_gblock_end-largest_gblock_start)
                largest_gblock_start = gblock_start
                largest_gblock_end = gblock_end
            end
        end

        # Write progress to the screen
        if (current_gpt % 1000 == 0) 
            @printf("Number of roots found up to Gram point %d: %d.\n", current_gpt, tot_num_roots)
        end
    end

    # Write results to the output files

    # The total number of roots
    @printf(info_file,"The total number of roots between Gram points N = %d and N = %d is %d.\n",
                            START,current_gpt,tot_num_roots)
    
    # The total number of Gram blocks
    @printf(info_file,"The total number of Gram blocks found is %d.\n", tot_num_gblocks)

    # The total number of roots in the Gram blocks
    @printf(info_file,"The total number of roots in the %d Gram blocks is %d.\n", 
                            tot_num_gblocks,tot_num_gblock_roots)

    # The longest Gram block found
    @printf(info_file,"The longest Gram block found is [%d, %d].\n",largest_gblock_start,
                            largest_gblock_end)

    # Close the file stream
    close(info_file)
end

main()