r"""
Reconstruction of algebraic forms from invariants

This module lists classical invariants and covariants of homogeneous
polynomials (also called algebraic forms) under the action of the
special linear group. That is, we are dealing with polynomials of
degree `d` in `n` variables. The special linear group `SL(n,\CC)` acts
on the variables `(x_1,\dots, x_n)` linearly,

.. MATH::

    (x_1,\dots, x_n)^t \to A (x_1,\dots, x_n)^t
    ,\qquad
    A \in SL(n,\CC)

The linear action on the variables transforms a polynomial `p`
generally into a different polynomial `gp`. We can think of it as an
action on the space of coefficients in `p`. An invariant is a
polynomial in the coefficients that is invariant under this action. A
covariant is a polynomial in the coefficients and the variables
`(x_1,\dots, x_n)` that is invariant under the combined action.

For example, the binary quadratic `p(x,y) = a x^2 + b x y + c y^2`
has as its invariant the discriminant `\mathop{disc}(p) = b^2 - 4 a
c`. This means that for any `SL(2,\CC)` coordinate change

.. MATH::

    \begin{pmatrix} x' \\ y' \end{pmatrix}
    =
    \begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}
    \begin{pmatrix} x \\ y \end{pmatrix}
    \qquad
    \alpha\delta-\beta\gamma=1

the discriminant is invariant, `\mathop{disc}\big(p(x',y')\big) =
\mathop{disc}\big(p(x,y)\big)`.

To use this module, you should use the factory object
:class:`invariant_theory <InvariantTheoryFactory>`. For example, take
the quartic::

    sage: R.<x,y> = QQ[]
    sage: q = x^4 + y^4
    sage: quartic = invariant_theory.binary_quartic(q);  quartic
    Binary quartic with coefficients (1, 0, 0, 0, 1)


One invariant of a quartic is known as the Eisenstein
D-invariant. Since it is an invariant, it is a polynomial in the
coefficients (which are integers in this example)::

    sage: quartic.EisensteinD()
    1

One example of a covariant of a quartic is the so-called g-covariant
(actually, the Hessian). As with all covariants, it is a polynomial in
`x`, `y` and the coefficients::

    sage: quartic.g_covariant()
    -x^2*y^2

As usual, use tab completion and the online help to discover the
implemented invariants and covariants.

In general, the variables of the defining polynomial cannot be
guessed. For example, the zero polynomial can be thought of as a
homogeneous polynomial of any degree. Also, since we also want to
allow polynomial coefficients we cannot just take all variables of the
polynomial ring as the variables of the form. This is why you will
have to specify the variables explicitly if there is any potential
ambiguity. For example::

    sage: invariant_theory.binary_quartic(R.zero(), [x,y])
    Binary quartic with coefficients (0, 0, 0, 0, 0)

    sage: invariant_theory.binary_quartic(x^4, [x,y])
    Binary quartic with coefficients (0, 0, 0, 0, 1)

    sage: R.<x,y,t> = QQ[]
    sage: invariant_theory.binary_quartic(x^4 + y^4 + t*x^2*y^2, [x,y])
    Binary quartic with coefficients (1, 0, t, 0, 1)

Finally, it is often convenient to use inhomogeneous polynomials where
it is understood that one wants to homogenize them. This is also
supported, just define the form with an inhomogeneous polynomial and
specify one less variable::

    sage: R.<x,t> = QQ[]
    sage: invariant_theory.binary_quartic(x^4 + 1 + t*x^2, [x])
    Binary quartic with coefficients (1, 0, t, 0, 1)

REFERENCES:

.. [WpInvariantTheory] :wikipedia:`Glossary_of_invariant_theory`
.. _[Cle1872] 

AUTHORS:

- Jesper Noordsij (2018-05-24): initial version

"""

#*****************************************************************************
#     Copyright (C) 2018 Jesper Noordsij <jesper.noordsij@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************






######################################################################

def binary_form_from_invariants(degree, invariants, type='default'):
    """
    Function to reconstruct a binary form from the values of its
    invariants.

    INPUT:

    - ``degree`` --  The degree of the binary form.

    - ``invariants`` --  The values of the invariants of the binary form.

    - ``type`` -- The type of invariants given.

    OUTPUT:

    A set of coefficients of a binary form, whose invariants are equal
    to the given ``invariants`` up to a scaling.

    EXAMPLES::

        sage: invariants = [1, 0, 0]
        sage: binary_form_from_invariants(5, invariants)
        (1, 0, 0, 0, 0, 1)

        sage: binary_form_from_invariants(6, invariants)
        Traceback (most recent call last):
        ...
        NotImplementedError: No reconstruction for binary forms of degree 6 implemented.


    """
    if degree == 5:
        return binary_quintic_from_invariants(invariants, type)
    else:
        raise NotImplementedError('No reconstruction for binary forms of degree %s implemented.' % degree)


######################################################################

def binary_quadratic_from_invariants(discriminant):
    """
    Function to reconstruct a binary quadratic from the value of its
    discriminant.

    """
    if discriminant == 0:
        return (1, 0, 0)
    else:
        return (1, 1, 0)


######################################################################

def binary_cubic_from_invariants(discriminant):
    """
    Function to reconstruct a binary cubic from the value of its
    discriminant.

    """
    if discriminant == 0:
        raise NotImplementedError('No distinction implemented for binary cubics with a double root.')
    else:
        return (1, 0, 0, 1)


######################################################################

def binary_quintic_from_invariants(invariants, type='clebsch', K=None, scaled=False, reduced=False):
    """
    Function to reconstruct a binary quintic from the values of its
    (Clebsch) invariants.

    INPUT:

    - ``invariants`` --  The values of the three or four invariants
      of the binary quintic.

    - ``type`` -- The type of invariants given. By default the
      given invariants are the invariants A, B, C (and R) as
      described by Clebsch in _[Cle1872].

    - ``K`` -- The field over which the quintic is defined.

    - ``scaled`` -- A boolean to determine wether the coefficients should
      be scaled so the result is independent of the scaling of the invariants.

    OUTPUT:

    A set of coefficients of a binary quintic, whose Clebsch invariants
    are equal to ``invariants`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_quintic_from_invariants
        sage: R.<x0, x1> = QQ[]
        sage: p = 3*x1^5 + 6*x1^4*x0 + 3*x1^3*x0^2 + 4*x1^2*x0^3 - 5*x1*x0^4 + 4*x0^5
        sage: quintic = invariant_theory.binary_quintic(p, x0, x1)
        sage: invs = [quintic.A_invariant(),quintic.B_invariant(),quintic.C_invariant()]
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
        (-9592267437341790539005557/244140625000000,
         -2149296928207625556323004064707/610351562500000000,
         -11149651890347700974453304786783/76293945312500000,
         -122650775751894638395648891202734239/47683715820312500000,
         -323996630945706528474286334593218447/11920928955078125000,
         -1504506503644608395841632538558481466127/14901161193847656250000)

    The form obtained corresponds with the one found when setting the coordinates
    equal to the covariants alpha and beta::

        sage: alpha = quintic.alpha_covariant()
        sage: beta = quintic.beta_covariant()
        sage: g = matrix([[alpha(x0=1,x1=0),alpha(x0=0,x1=1)],[beta(x0=1,x1=0),beta(x0=0,x1=1)]])^-1
        sage: g = g*g.determinant()^-1
        sage: quintic.transformed(g).coeffs() == reconstructed
        True

    We can check that the invariants match by scaling the invariants
    B and C to see if they match::

        sage: newform = sum([ reconstructed[i]*x0^i*x1^(5-i) for i in range(6) ])
        sage: newquintic = invariant_theory.binary_quintic(newform, x0, x1)
        sage: scale = invs[0]/newquintic.A_invariant()
        sage: invs[1] == newquintic.B_invariant()*scale^2
        True
        sage: invs[2] == newquintic.C_invariant()*scale^3
        True

    If the invariant M vanishes, the coefficients are computed in a
    different way::

        sage: reconstructed = binary_quintic_from_invariants([3,1,2], K=QQ)
        sage: reconstructed
        (-66741943359375/2097152,
         125141143798828125/134217728,
         0,
         -52793920040130615234375/34359738368,
         19797720015048980712890625/1099511627776,
         4454487003386020660400390625/17592186044416)
        sage: newform = sum([ reconstructed[i]*x0^i*x1^(5-i) for i in range(6) ])
        sage: newquintic = invariant_theory.binary_quintic(newform, x0, x1)
        sage: scale = 3/newquintic.A_invariant()
        sage: [3, newquintic.B_invariant()*scale^2, newquintic.C_invariant()*scale^3]
        [3, 1, 2]

    Several special cases::

        sage: quintic = invariant_theory.binary_quintic(x0^5 - x1^5, x0, x1)
        sage: invs = [quintic.A_invariant(),quintic.B_invariant(),quintic.C_invariant()]
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
        (1, 0, 0, 0, 0, 1)
        sage: quintic = invariant_theory.binary_quintic(x0*x1*(x0^3-x1^3), x0, x1)
        sage: invs = [quintic.A_invariant(),quintic.B_invariant(),quintic.C_invariant()]
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
         (0, 1, 0, 0, 1, 0)

    For fields of characteristic 2, 3 or 5, there is no reconstruction implemented::

        sage: binary_quintic_from_invariants([3,1,2], K=GF(5))
        Traceback (most recent call last):
        ...
        NotImplementedError: No reconstruction implemented for fields of characteristic 2, 3 or 5.

    """
    if reduced:
        if len(invariants) == 3:
            invariants = reduce_invariants(invariants, [1,2,3])
        else:
            invariants = reduce_invariants(invariants, [2,4,6,9])
    A, B, C = invariants[0:3]
    if K is None:
        K = A.parent()
    if K.characteristic() in [2, 3, 5]:
        raise NotImplementedError('No reconstruction implemented for fields of characteristic 2, 3 or 5.')
    M = 2*A*B - 3*C
    N = K(2)**-1 * (A*C-B**2)
    R2 = -K(2)**-1 * (A*N**2-2*B*M*N+C*M**2)
    scale = [1,1,1,1,1,1]
    from sage.functions.all import binomial, sqrt
    try:
        if R2.is_square():
            R = sqrt(R2)
        else:
            # if R2 is not a square, we scale the invariants in a suitable way so that the 'new' R2 is a square
            # r = R2.squarefree_part() # slow!
            invariants = [R2*A, R2**2*B, R2**3*C]
            return binary_quintic_from_invariants(invariants, type, K, scaled, reduced)
    except (AttributeError, NotImplementedError):
        if len(invariants) > 3:
            R = invariants[3]
        else:
            raise ValueError('Value of R could not be determined.')
    if M == 0:
        if N == 0:
            if A == 0:
                raise NotImplementedError('No reconstruction implemented for quintics with a treefold linear factor.')
            else:
                if B == 0:
                    return (1,0,0,0,0,1) # x**5 + y**5
                else:
                    return (0,1,0,0,1,0) # x*y*(x**3+y**3)
        else:
            # case corresponding to using alpha and gamma as coordinates
            if A == 0:
                return (1,0,0,0,1,0) # x*(x**4+y**4)
            else:
                if scaled:
                    scale = [ A**-14*(R/A**3)**i for i in range(6) ] # subs = {y:(R/A**3)*y}
                D = -N
                Delta = C
                A1 = (2*K(3)**-1*A**2-B)*N*B*K(2)**-1 - N**2*K(2)**-1
                B0 = 2*K(3)**-1*A*R
                B1 = A*N*B*K(3)**-1
                C0 = 2*K(3)**-1*R
                C1 = B*N
    else:
        # case corresponding to using alpha and beta as coordinates
        if R == 0:
            if A == 0:
                return (1,0,10,0,-15,0) # x**5 + 10*x**3*y**2 - 15*x*y**4
            elif scaled:
                scale = [ sqrt(A)**(i-18) for i in range(6) ] # subs = {x:x/A, y:y/sqrt(A)}
        else:
            if A == 0:
                if B == 0:
                    return (1,0,0,1,0,0) # x**2*(x**3+y**3)
                elif scaled:
                    scale = [ R**-2*(R/B**2)**i for i in range(6) ] # subs = {y:(R/B**2)*y}
            elif scaled:
                scale = [ A**-9*(R/A**4)**i for i in range(6) ] # subs = {y:(R/A**4)*y}
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
    # D**(-5)*(A5*y**5 - 5*A4*x*y**4 + 10*A3*x**2*y**3 - 10*A2*x**3*y**2 + 5*A1*x**4*y- A0*x**5)
    #tuple([K((-1)**(i+1)*D**(-5)*binomial(5,i)*scale[5-i]*eval('A%d' %i)) for i in range(6)])
    coeffs = tuple([K((-1)**(i+1)*binomial(5,i)*scale[5-i]*eval('A%d' %i)) for i in range(6)])
    if reduced:     
        from sage.arith.misc import gcd
        return tuple([coeffs[i]/gcd(coeffs) for i in range(6)])
    else:
        return coeffs

        
def reduce_invariants(invariants, weights):
    factors = [dict(I.factor()) for I in invariants]
    scalar = 1
    n = len(weights)
    from sage.arith.misc import gcd
    for prime in gcd(invariants).factor():
        p = prime[0]
        for D in factors:
            if not D.has_key(p):
                D[p] = 0
        scalar = scalar*p**min([factors[i][p]//weights[i] for i in range(n)])
        print scalar.parent(); print scalar
    return [invariants[i]*scalar**-weights[i] for i in range(n)]

