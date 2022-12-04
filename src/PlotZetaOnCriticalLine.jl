#=**
**
**    Plot (t,Z(t)) through the Riemann-Siegel formula.
**
**************************************************************************
**    Matthew Kehoe
**    11/06/2022
**
**    This program uses the Riemann-Siegel formula to compute approximations
**    to Z(t). The program returns coordinate pairs (t,Z(t)) for a <= t <= b 
**    and plots the results through the standard plotting tools in Julia.
**
**    Calculations should be performed for t > 10.
**************************************************************************=#

using Plots
using Figures
using LaTeXStrings
using Printf

# Use the functions defined in zetaFunctions.jl
include("zetaFunctions.jl")

function main()

    a = 8            # lower bound of t domain
    b = 30.0               # upper bound of t domain
    numSamples = 500         # number of samples to compute

    samples = numSamples - 1
    sampleIndex = 0 
    intervalLength = b - a

    rePart = Float64[]
    imPart = Float64[]
    axisVal = Float64[]

    while sampleIndex <= samples
        t = a + sampleIndex/samples * intervalLength
        reZeta = RiemennZ(t,4)*cos(theta(t))
        imZeta = -1.0*RiemennZ(t,4)*sin(theta(t))

        @printf("t = %.12f\t zeta(1/2+it) = %.12f + i %.12f\n",t,reZeta,imZeta);

        push!(rePart, reZeta)
        push!(imPart, imZeta)
        push!(axisVal, t)

        sampleIndex = sampleIndex + 1
    end

    # The plot below is the zeta function spiral
    display(plot((rePart, imPart), color=:blue, width=3, title=L"\zeta(1/2 + it)", 
          xaxis= L"Re(\zeta(1/2 + it))", yaxis= L"Im(\zeta(1/2 + it))", 
          framestyle= :origin,label="RiemannZ"))
    #scatter!((rePart, imPart),color=:red, label="Values")

    plot((axisVal,rePart), color=:blue, width=3, title="Real and Imaginary Parts on Critical Line", xaxis=L"s", 
          yaxis = L"\zeta(1/2 + it)", framestyle= :origin,label=L"Re(\zeta(s))")

    plot!((axisVal,imPart), color=:red, width=3, title="Real and Imaginary Parts on Critical Line", xaxis=L"s", 
          yaxis = L"\zeta(1/2 + it)", framestyle= :origin,label=L"Im(\zeta(s))")
end

main()