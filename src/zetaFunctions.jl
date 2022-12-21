#=**
**
**    Assorted functions for use with other Zeta programs.
**
**************************************************************************
**    Matthew Kehoe
**    10/30/2022
**
**    Assortment of different functions for use with other Zeta programs.
**    Contains both the Riemann-Siegel Z(t) function alongside the Riemann-
**    Siegel theta function and its derivative.
**************************************************************************=#

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
Returns -1 if the argument is odd and +1 if the argument is even.
### Input
@param n - an integer.
### Output
@return -1 if the integer is odd and +1 if the integer is even
"""
function even(n::Int)
    if n % 2 == 0
        return 1
    else 
        return -1
    end
end


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

  while abs(tn - tn1) > 1e-8
      tn = tn1
      # Compute t_1 = t_0 - f(t_0) / f'(t_0) 
      # Here, f is the Riemann-Siegel theta function minus n times pi
      tn1 = tn - (theta(tn) - n*pi) / thetaDeriv(tn)
  end

  return tn1
end


"""
    Riemann-Siegel Z(t) function implemented per the Riemenn Siegel
    formula. See http://mathworld.wolfram.com/Riemann-SiegelFormula.html
    for details.
### Input
    @param t - the value of t inside the Z(t) function.
    @param r - referenced for calculating the remainder terms by the
         Taylor series approximations.
### Output
    @return the approximate value of Z(t) through the Riemann-Siegel
         formula
"""
function RiemennZ(t::Float64, r::Int)
    tDivTwoPi = t / (2.0 * pi)
    val = sqrt(tDivTwoPi)
    N = floor(abs(val))
    sum = 0.0

    i = 1
    while i <= N
          sum += (cos(theta(t) - t * log(i))) / sqrt(i)
          i += 1
    end
    sum = 2.0 * sum

    frac = val - N
    k = 0
    R = 0.0

    # Calculate every remainder term through the Taylor Series coefficients.
    # The coefficients are defined in the Riemann-Siegel C function.
    while k <= r
        R += C(k, 2.0*frac-1.0) * tDivTwoPi^(k * (-0.5))
        k += 1
    end

    remainder = (-1)^(N-1) * tDivTwoPi^(-0.25) * R
    return sum + remainder
end


"""
    C terms for the Riemann-Siegel formula. See
    Haselgrove (1960): Tables of the Riemann Zeta Function,
    doi:10.2307/3614148  for details.
    Calculates the Taylor Series coefficients for C0, C1, C2, C3,
    and C4.
### Input
    @param n - the number of coefficient terms to use.
    @param z - referenced per the Taylor series calculations.
### Output
    @return the Taylor series approximation of the remainder terms.
"""
function C(n::Int, z::Float64)
    if n==0
        return (.38268343236508977173 * z^0.0
            +.43724046807752044936 * z^2.0
            +.13237657548034352332 * z^4.0
            -.01360502604767418865 * z^6.0
            -.01356762197010358089 * z^8.0
            -.00162372532314446528 * z^10.0
            +.00029705353733379691 * z^12.0
            +.00007943300879521470 * z^14.0
            +.00000046556124614505 * z^16.0
            -.00000143272516309551 * z^18.0
            -.00000010354847112313 * z^20.0
            +.00000001235792708386 * z^22.0
            +.00000000178810838580 * z^24.0
            -.00000000003391414390 * z^26.0
            -.00000000001632663390 * z^28.0
            -.00000000000037851093 * z^30.0
            +.00000000000009327423 * z^32.0
            +.00000000000000522184 * z^34.0
            -.00000000000000033507 * z^36.0
            -.00000000000000003412 * z^38.0
            +.00000000000000000058 * z^40.0
            +.00000000000000000015 * z^42.0)
    elseif n==1
        return (-.02682510262837534703 * z^1.0
            +.01378477342635185305 * z^3.0
            +.03849125048223508223 * z^5.0
            +.00987106629906207647 * z^7.0
            -.00331075976085840433 * z^9.0
            -.00146478085779541508 * z^11.0
            -.00001320794062487696 * z^13.0
            +.00005922748701847141 * z^15.0
            +.00000598024258537345 * z^17.0
            -.00000096413224561698 * z^19.0
            -.00000018334733722714 * z^21.0
            +.00000000446708756272 * z^23.0
            +.00000000270963508218 * z^25.0
            +.00000000007785288654 * z^27.0
            -.00000000002343762601 * z^29.0
            -.00000000000158301728 * z^31.0
            +.00000000000012119942 * z^33.0
            +.00000000000001458378 * z^35.0
            -.00000000000000028786 * z^37.0
            -.00000000000000008663 * z^39.0
            -.00000000000000000084 * z^41.0
            +.00000000000000000036 * z^43.0
            +.00000000000000000001 * z^45.0)
    elseif n==2
        return(+.00518854283029316849 * z^0.0
            +.00030946583880634746 * z^2.0
            -.01133594107822937338 * z^4.0
            +.00223304574195814477 * z^6.0
            +.00519663740886233021 * z^8.0
            +.00034399144076208337 * z^10.0
            -.00059106484274705828 * z^12.0
            -.00010229972547935857 * z^14.0
            +.00002088839221699276 * z^16.0
            +.00000592766549309654 * z^18.0
            -.00000016423838362436 * z^20.0
            -.00000015161199700941 * z^22.0
            -.00000000590780369821 * z^24.0
            +.00000000209115148595 * z^26.0
            +.00000000017815649583 * z^28.0
            -.00000000001616407246 * z^30.0
            -.00000000000238069625 * z^32.0
            +.00000000000005398265 * z^34.0
            +.00000000000001975014 * z^36.0
            +.00000000000000023333 * z^38.0
            -.00000000000000011188 * z^40.0
            -.00000000000000000416 * z^42.0
            +.00000000000000000044 * z^44.0
            +.00000000000000000003 * z^46.0)
    elseif n==3
        return (-.00133971609071945690 * z^1.0
            +.00374421513637939370 * z^3.0
            -.00133031789193214681 * z^5.0
            -.00226546607654717871 * z^7.0
            +.00095484999985067304 * z^9.0
            +.00060100384589636039 * z^11.0
            -.00010128858286776622 * z^13.0
            -.00006865733449299826 * z^15.0
            +.00000059853667915386 * z^17.0
            +.00000333165985123995 * z^19.0
            +.00000021919289102435 * z^21.0
            -.00000007890884245681 * z^23.0
            -.00000000941468508130 * z^25.0
            +.00000000095701162109 * z^27.0
            +.00000000018763137453 * z^29.0
            -.00000000000443783768 * z^31.0
            -.00000000000224267385 * z^33.0
            -.00000000000003627687 * z^35.0
            +.00000000000001763981 * z^37.0
            +.00000000000000079608 * z^39.0
            -.00000000000000009420 * z^41.0
            -.00000000000000000713 * z^43.0
            +.00000000000000000033 * z^45.0
            +.00000000000000000004 * z^47.0)
    else
        return (+.00046483389361763382 * z^0.0
            -.00100566073653404708 * z^2.0
            +.00024044856573725793 * z^4.0
            +.00102830861497023219 * z^6.0
            -.00076578610717556442 * z^8.0
            -.00020365286803084818 * z^10.0
            +.00023212290491068728 * z^12.0
            +.00003260214424386520 * z^14.0
            -.00002557906251794953 * z^16.0
            -.00000410746443891574 * z^18.0
            +.00000117811136403713 * z^20.0
            +.00000024456561422485 * z^22.0
            -.00000002391582476734 * z^24.0
            -.00000000750521420704 * z^26.0
            +.00000000013312279416 * z^28.0
            +.00000000013440626754 * z^30.0
            +.00000000000351377004 * z^32.0
            -.00000000000151915445 * z^34.0
            -.00000000000008915418 * z^36.0
            +.00000000000001119589 * z^38.0
            +.00000000000000105160 * z^40.0
            -.00000000000000005179 * z^42.0
            -.00000000000000000807 * z^44.0
            +.00000000000000000011 * z^46.0
            +.00000000000000000004 * z^48.0)
    end
end
