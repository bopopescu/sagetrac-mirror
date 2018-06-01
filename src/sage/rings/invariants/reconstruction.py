r"""
Reconstruction of algebraic forms from invariants



AUTHORS:

- Jesper Noordsij (2018-05-24): initial version

"""

#*****************************************************************************
#       Copyright (C) 2018 Jesper Noordsij <jesper.noordsij@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************






######################################################################

def binary_form_from_invariants(degree, invariants, type='default', as_form=False):
    if degree == 5:
        coefficients = binary_quintic_from_invariants(invariants, type)
        if as_form == True:
            return invariant_theory.binary_quintic(as_form) # TODO: fix variables!@!
        else:
            return coefficients
    else:
        return NotImplementedError('No reconstruction for binary forms of degree %s implemented.' % degree)
    
    
######################################################################

def binary_quintic_from_invariants(invariants, type, K=None):
    r"""
    Function to reconstruct a binary quintic from the values of its 
    Clebsch invariants

    INPUT:

    - ``invs`` --  The values of the three Clebsch invariant A, B
      and C of the binary quintic.

    - ``x, y`` -- The variables for the resulting binary quintic.
    
    - ``K`` -- The field over which the quintic will be defined.

    OUTPUT:

    A binary quintic whose Clebsch invariants are equal to ``invs``,
    up to a scaling of the invariants.

    EXAMPLES::
        
        sage: K = GF(7)
        sage: R.<x0,x1> = K[]
        sage: p = 3*x1^5 + 6*x1^4*x0 + 3*x1^3*x0^2 + 4*x1^2*x0^3 - 5*x1*x0^4 + 4*x0^5
        sage: quintic = invariant_theory.binary_quintic(p, x0, x1)
        sage: invs = [K(quintic.A_invariant()),K(quintic.B_invariant()),K(quintic.C_invariant())]
        sage: reconstructed = quintic_from_clebsch_invariants(invs, x0, x1)
        sage: reconstructed
        -2*x0^5 + 3*x0^4*x1 + x0^3*x1^2 + 3*x0^2*x1^3 + x0*x1^4 + 3*x1^5
    
    The form obtained corresponds with the one found when setting the coordinates
    equal to the covariants alpha and beta::
                
        sage: alpha = quintic.alpha_covariant()
        sage: beta = quintic.beta_covariant()
        sage: g = matrix([[alpha(x0=1,x1=0),alpha(x0=0,x1=1)],[beta(x0=1,x1=0),beta(x0=0,x1=1)]])^-1
        sage: quintic.transformed(g).form() == reconstructed
        True
        
    We can check that the invariants match by scaling the invariants
    B and C to see if they match::
    
        sage: newquintic = invariant_theory.binary_quintic(reconstructed, x0, x1)
        sage: scale = invs[0]/newquintic.A_invariant()
        sage: invs[1] == newquintic.B_invariant()*scale^2
        True
        sage: invs[2] == newquintic.C_invariant()*scale^3
        True
    
    If the invariant M vanishes, the coefficients are computed in a slightly
    different way::
    
        sage: reconstructed = quintic_from_clebsch_invariants([K(3),1,2], x0, x1)
        sage: reconstructed
        x0^5 + 2*x0^4*x1 - 2*x0^3*x1^2 - x0*x1^4 - x1^5
        sage: newquintic = invariant_theory.binary_quintic(reconstructed, x0, x1)
        sage: scale = 3/newquintic.A_invariant()
        sage: [3, newquintic.B_invariant()*scale^2, newquintic.C_invariant()*scale^3]
        [3, 1, 2]
        
    Several special cases::
    
        sage: quintic = invariant_theory.binary_quintic(x0^5 - x1^5, x0, x1)
        sage: invs = [K(quintic.A_invariant()),K(quintic.B_invariant()),K(quintic.C_invariant())]
        sage: reconstructed = quintic_from_clebsch_invariants(invs, x0, x1)
        sage: reconstructed
        x0^5 - x1^5
        sage: quintic = invariant_theory.binary_quintic(x0*x1*(x0^3-x1^3), x0, x1)
        sage: invs = [K(quintic.A_invariant()),K(quintic.B_invariant()),K(quintic.C_invariant())]
        sage: reconstructed = quintic_from_clebsch_invariants(invs, x0, x1)
        sage: reconstructed
        x0^4*x1 - x0*x1^4

    """
    A, B, C = invariants
    if K is None:
        K = A.parent()
    if K.characteristic() in [2, 3, 5]:
        return NotImplementedError('No reconstruction implemented for fields of characteristic 2, 3 or 5.')
    M = 2*A*B - 3*C
    N = K(2)**-1 * (A*C-B**2)
    R2 = -K(2)**-1 * (A*N**2-2*B*M*N+C*M**2)
    from sage.functions.all import binomial, sqrt
    if R2.is_square():
        R = sqrt(R2)
    else:
        # if R2 is not a square, we scale the invariants in a suitable way so that the 'new' R2 is a square
        invariants = [R2*A, R2**2*B, R2**3*C]
        return binary_quintic_from_invariants(invariants, type, K)
    if M == 0:
        if N == 0:
            if A == 0:
                return NotImplementedError('No reconstruction implemented for quintics with a treefold linear factor.')
            else:
                if B == 0:
                    return (1,0,0,0,0,1) # x**5 + y**5
                else:
                    return (0,1,0,0,1,0) # x*y*(x**3+y**3)
        else:
            # case corresponding to using alpha and gamma as coordinates
            D = -N
            Delta = C
            A1 = (2*K(3)**-1*A**2-B)*N*B*K(2)**-1 - N**2*K(2)**-1
            B0 = 2*K(3)**-1*A*R
            B1 = A*N*B*K(3)**-1
            C0 = 2*K(3)**-1*R
            C1 = B*N
    else:
        # case corresponding to using alpha and beta as coordinates
        D = -M
        Delta = A
        A1 = (2*K(3)**-1*A**2-B)*(N*A-M*B)*K(2)**-1 - M*(N*K(2)**-1-M*A*K(3)**-1)
        B0 = R
        B1 = K(2)**-1*(N*A-M*B)
        C0 = 0
        C1 = -M
    A0 = (2*K(3)**-1*A**2-B)*R
    A2 = -D*B0 - K(2)**-1*Delta*A0
    A3 = -D*B1 - K(2)**-1*Delta*A1
    A4 = D**2*C0 + D*Delta*B0 + K(4)**-1*Delta**2*A0
    A5 = D**2*C1 + D*Delta*B1 + K(4)**-1*Delta**2*A1
    return tuple([D**(-5)*binomial(6,i)*eval('A%d' %i) for i in range(6)]) # D**(-5)*(A5*x**5 + 5*A4*x**4*y + 10*A3*x**3*y**2 + 10*A2*x**2*y**3 + 5*A1*x*y**4 + A0*y**5)