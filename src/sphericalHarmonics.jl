"""
https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
"""

using Base.Test
using Polynomials

"""
Calculate the Legendre polynomial P_n as a Poly polynomial.
"""
function SH_Pn( n::Int )
    return sum([ binomial(n,k)*binomial(n+k,k)*Poly([-1/2,1/2])^k for k in 0:n ])
end  #function

@test SH_Pn( 3 ) == Poly( [0,-1.5,0,2.5] )
@test SH_Pn( 0 ) == Poly( [1] )
@test SH_Pn( -1 ) == Poly( [0] )

"""
Calculate the value of Legendre polynomial P_n at x.
"""
function SH_Pn( n::Int,x::Number )
	return polyval( SH_Pn(n),x )
end  #function

@test SH_Pn( 3,-1 ) == -1.0
@test SH_Pn( 3,0 )  ==  0.0
@test SH_Pn( 3,+1 ) == +1.0

"""
Calculate the value of the associated Legendre polynomial P_l^m at x.

Guaranteed to work only in [-1,+1].
"""
function SH_Plm( l::Int,m::Int,x )
	return (-1)^m*(1-x^2)^(m/2)*polyder( SH_Pn(l),m )(x)
end  #function

@test SH_Plm( 1,1,-1 )   ≈  0.0
@test SH_Plm( 1,1,0 )    ≈ -1.0
@test SH_Plm( 1,1,+1 )   ≈  0.0
@test SH_Plm( 2,2,-1/2 ) ≈  2.25
@test SH_Plm( 2,2,0 )    ≈  3.0
@test SH_Plm( 2,2,+1/2 ) ≈  2.8125
@test SH_Plm( 3,1,-1/2 ) ≈ -0.3247595264191645
@test SH_Plm( 3,1,0 )    ≈  1.5
@test SH_Plm( 3,1,+1/2 ) ≈ -0.3247595264191645

"""
Calculate the value of spherical harmonics Y_l^m(θ,φ).

The quantum-mechanical normalization convention is used.
"""
function Ylm(l,m,θ,φ)
    N = (-1)^m*√((2*l+1)/(4π)*factorial(l-m)/factorial(l+m))
	return N*exp(1im*φ)*SH_Plm(l,m,cos(θ))
end  #function

@test Ylm( 0,0,0,0 )     ≈  1/2*√(1/pi)
@test Ylm( 1,1,π/2,0 )   ≈  1/2*√(3/(2π))
@test Ylm( 1,0,0,0 )     ≈  1/2*√(3/π)
@test Ylm( 3,1,π/2,0 )   ≈ -1/8*√(21/π)
@test Ylm( 4,2,π/4,π/4 ) ≈  3/8*√(5/(2π))*√2/4*(1+1im)*(7/2-1)

