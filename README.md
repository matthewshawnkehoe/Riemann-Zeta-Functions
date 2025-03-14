# Riemann Siegel and Assorted Zeta Functions
The German mathematician Bernhard Riemann published his work on the investigation of prime numbers in 1859. This eight-page paper, titled "On the Number of Primes Less Than a Given Magnitude", expanded on the work of Euler, Gauss, and Dirichlet and made numerous contributions to the study of both complex and analytic function theory. One of the major claims of this work was a conjecture about the location of zeros of the Riemann zeta function. To put it simply, Riemann hypothesized that all non-trivial zeros of the Riemann zeta function must have real part equal to $1/2$ and lie inside the critical strip $0 < s < 1$. This conjecure is now known as the Riemann Hypothesis (RH). Numerous mathematicians and computer scientists have created numerical algorithms to verify the Riemann Hypothesis up to a specific height $T$. The three most famous of these numerical algorithms are:

1. Euler-Maclaurin Summation Formula
2. Riemann–Siegel Formula
3. Odlyzko–Schonhage Algorithm

This repository contains a Julia implementation of the Riemann–Siegel formula. The code was written by Matthew Kehoe (mskehoe001@gmail.com) and is partially based off the work (written in C) by [Glen Pugh](https://web.viu.ca/pughg/) and [Ken Takusagawa](http://web.mit.edu/kenta/www/six/parallel/2-Final-Report.html). Future commits of this repository may contain a working implementation of the Odlyzko–Schonhage algorithm in Julia and various enhancements to the root-finding and plotting programs.

## Installation
Download the Julia code in the src directory. The individual programs are

1. `BernoulliNumbers.jl:` This program applies the Akiyama-Tanigawa algorithm to display integer values of the Riemann zeta function.
2. `PlotZetaOffCriticalLine.jl:` Calculates real values of the Riemann zeta function off the critical line and plots them to the screen.
3. `PlotZetaOnCriticalLine.jl:` Uses the Riemann-Siegel formula to plot the Riemann-Siegel $Z$ function, $(t,Z(t))$, at specified time intervals.
4. `RiemannSiegel.jl:` Implementation of the Riemann-Siegel formula to calculate non-trivial zeros of the Riemann zeta function.
5. `TuringTest.jl:` Applies Turing's method to bound $S(g_n)$ in the expression $(-1)^n  Z(g_n + h_n) > 0$ (which can be used to verify the RH up to a specific height $T$).
6. `rootCount.jl:` Verifies the RH up to a final Gram point using the Riemann-Siegel formula.
7. `zetaFunctions.jl:` A library of functions used in all the other programs.

## Plotting 

Example 1: Plot of $\zeta(1/2+it)$ on the critical line.

![image](https://github.com/matthewshawnkehoe/Riemann-Zeta-Functions/assets/13266088/9c7cdfbe-135e-45ef-822f-65c91723ced3)

Example 2: Real and Imaginary values of the Riemann-Siegel $Z$ function.

![image](https://github.com/matthewshawnkehoe/Riemann-Zeta-Functions/assets/13266088/9629886e-2eed-4f34-ba3a-c35f8919c418)

Example 3: Real values of the Riemann Zeta function.

![image](https://github.com/matthewshawnkehoe/Riemann-Zeta-Functions/assets/13266088/de59d50f-35f5-475e-bfac-f05330a10383)


## References
Interested readers can review the following material:

* Overview of the Riemann Zeta Function: [Edward's Riemann Zeta Function](https://www.amazon.com/Riemanns-Zeta-Function-Harold-Edwards/dp/0486417409)
* C Implementation of Riemann-Siegel Formula: [Glen's Project](https://web.viu.ca/pughg/thesis.d/masters.thesis.pdf) and [Ken's Project](http://web.mit.edu/kenta/www/six/parallel/2-Final-Report.html)
* Beamer Presentation: [Calculating zeros of the RZF]([https://axion004.files.wordpress.com/2022/12/calculating_zeros_of_the_riemann_zeta_function.pdf](https://matthewshawnkehoe.github.io/files/kehoe_calculating_zeros_of_the_riemann_zeta_function.pdf)

## Future Work

1. Implement the Odlyzko–Schonhage Algorithm in Julia.
2. Apply more remainder terms in the calculation of the Riemann-Siegel formula (Ken computed a lot more than four).
3. Implement 3D plots of the Riemann zeta function with real and complex values.

## Questions

Email Matthew Kehoe (mskehoe001@gmail.com) about questions related to this work.

