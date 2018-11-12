r"""
Constructions for common cones

This module provides a ``cones`` object with shortcut methods to create
a few types of cones:

  - The nonnegative orthant.

  - The rearrangement cone of order ``p``.

  - The Schur cone.

  - The trivial cone.

At the moment, only convex rational polyhedral cones are
supported---specifically, those cones that can be built using the
:func:`Cone` constructor. As a result, each shortcut method can be
passed either an ambient dimension ``ambient_dim``, or a toric
``lattice`` (from which the dimension can be inferred) to determine
the ambient space.

Here are some typical usage examples::

    sage: cones.nonnegative_orthant(2).rays()
    N(1, 0),
    N(0, 1)
    in 2-d lattice N

::

    sage: cones.rearrangement(2,2).rays()
    N( 1,  0),
    N( 1, -1),
    N(-1,  1)
    in 2-d lattice N

::

    sage: cones.schur(3).rays()
    N(1, -1,  0),
    N(0,  1, -1)
    in 3-d lattice N

::

    sage: cones.trivial(3).rays()
    Empty collection
    in 3-d lattice N

To specify some other lattice, pass it as an argument to the method::

    sage: K = cones.nonnegative_orthant(3)
    sage: cones.schur(lattice=K.dual().lattice())
    2-d cone in 3-d lattice M

For more information about these cones, see the documentation for the
individual methods on the :class:`ConvexRationalPolyhedralConeFactory`
class and the references therein.

"""

from __future__ import print_function
from six.moves import range

from cone import Cone
from sage.geometry.toric_lattice import ToricLattice
from sage.matrix.all import matrix
from sage.rings.all import ZZ



class ConvexRationalPolyhedralConeFactory:
    r"""
    Convenient shortcuts for common cones.

    .. WARNING::

        You should not instantiate this class. Use the
        existing ``cones`` object instead.

    """

    def nonnegative_orthant(self, ambient_dim=None, lattice=None):
        r"""
        The nonnegative orthant in ``ambient_dim`` dimensions, or living
        in ``lattice``.

        INPUT:

        - ``ambient_dim`` -- (default: ``None``) the dimension of the
          ambient space, a nonnegative integer.

        - ``lattice`` -- (default: ``None``) the toric lattice in which
          the cone will live.

        If ``ambient_dim`` is omitted, then it will be inferred from the
        rank of ``lattice``. If the ``lattice`` is omitted, then the
        default lattice of rank ``ambient_dim`` will be used. A
        ``ValueError`` is raised if neither ``ambient_dim`` nor
        ``lattice`` are specified.

        It is a ``ValueError`` to specify both ``ambient_dim`` and ``lattice``
        unless the rank of ``lattice`` is equal to ``ambient_dim``.

        OUTPUT:

        A :class:`.ConvexRationalPolyhedralCone` living in ``lattice``
        and having ``ambient_dim`` standard basis vectors as its
        generators. Each generating ray has the integer ring as its
        base ring.

        A ``ValueError`` is raised if,

        - Neither the ambient dimension nor the ambient lattice are specified.

        - Both the ambient dimension and the ambient lattice are specified,
          but they are incompatible.

        REFERENCES:

        - [BV2009]_

        EXAMPLES::

            sage: cones.nonnegative_orthant(3).rays()
            N(1, 0, 0),
            N(0, 1, 0),
            N(0, 0, 1)
            in 3-d lattice N

        TESTS:

        We can construct the trivial cone as the nonnegative orthant in a
        trivial vector space::

            sage: cones.nonnegative_orthant(0)
            0-d cone in 0-d lattice N

        The nonnegative orthant is a proper cone::

            sage: set_random_seed()
            sage: ambient_dim = ZZ.random_element(10)
            sage: K = cones.nonnegative_orthant(ambient_dim)
            sage: K.is_proper()
            True

        If a ``lattice`` was given, it is actually used::

            sage: L = ToricLattice(3, 'M')
            sage: cones.nonnegative_orthant(lattice=L)
            3-d cone in 3-d lattice M

        Unless the rank of the lattice disagrees with ``ambient_dim``::

            sage: L = ToricLattice(1, 'M')
            sage: cones.nonnegative_orthant(3, lattice=L)
            Traceback (most recent call last):
            ...
            ValueError: lattice rank=1 and ambient_dim=3 are incompatible

        We also get an error if no arguments are given::

            sage: cones.nonnegative_orthant()
            Traceback (most recent call last):
            ...
            ValueError: either the ambient dimension or the lattice must
            be specified

        """
        if ambient_dim is None and lattice is None:
            raise ValueError('either the ambient dimension or the lattice '
                             'must be specified')

        if ambient_dim is None:
            ambient_dim = lattice.rank()

        if lattice is None:
            lattice = ToricLattice(ambient_dim)

        if lattice.rank() != ambient_dim:
            raise ValueError('lattice rank=%d and ambient_dim=%d '
                             'are incompatible' % (lattice.rank(), ambient_dim))

        I = matrix.identity(ZZ,ambient_dim)
        return Cone(I.rows(), lattice)


    def rearrangement(self, p, ambient_dim=None, lattice=None):
        r"""
        The rearrangement cone of order ``p`` in ``ambient_dim``
        dimensions, or living in ``lattice``.

        The rearrangement cone of order ``p`` in ``ambient_dim``
        dimensions consists of all vectors of length ``ambient_dim``
        whose smallest ``p`` components sum to a nonnegative number.

        For example, the rearrangement cone of order one has its single
        smallest component nonnegative. This implies that all components
        are nonnegative, and that therefore the rearrangement cone of
        order one is the nonnegative orthant in its ambient space.

        When ``p`` and ``ambient_dim`` are equal, all components of the
        cone's elements must sum to a nonnegative number. In other
        words, the rearrangement cone of order ``ambient_dim`` is a
        half-space.

        INPUT:

        - ``p`` -- the number of components to "rearrange," an integer
          between ``1`` and ``ambient_dim`` inclusive.

        - ``ambient_dim`` -- (default: ``None``) the dimension of the
          ambient space, a nonnegative integer.

        - ``lattice`` -- (default: ``None``) the toric lattice in which
          the cone will live.

        If ``ambient_dim`` is omitted, then it will be inferred from the
        rank of ``lattice``. If the ``lattice`` is omitted, then the
        default lattice of rank ``ambient_dim`` will be used. A
        ``ValueError`` is raised if neither ``ambient_dim`` nor
        ``lattice`` are specified.

        It is a ``ValueError`` to specify both ``ambient_dim`` and ``lattice``
        unless the rank of ``lattice`` is equal to ``ambient_dim``.

        OUTPUT:

        A :class:`.ConvexRationalPolyhedralCone` representing the
        rearrangement cone of order ``p`` living in ``lattice``, with
        ambient dimension ``ambient_dim``. Each generating ray has the
        integer ring as its base ring.

        A ``ValueError`` is raised if,

        - Neither the ambient dimension nor the ambient lattice are specified.

        - Both the ambient dimension and the ambient lattice are specified,
          but they are incompatible.

        ALGORITHM:

        The generators for the rearrangement cone are given by [Jeong2017]_
        Theorem 5.2.3.

        REFERENCES:

        - [GJ2016]_

        - [HS2010]_

        - [Jeong2017]_

        EXAMPLES:

        The rearrangement cones of order one are nonnegative orthants::

            sage: orthant = cones.nonnegative_orthant(6)
            sage: cones.rearrangement(1,6).is_equivalent(orthant)
            True

        When ``p`` and ``ambient_dim`` are equal, the rearrangement cone
        is a half-space, so we expect its lineality to be one less than
        ``ambient_dim`` because it will contain a hyperplane but is not
        the entire space::

            sage: cones.rearrangement(5,5).lineality()
            4

        Jeong's Proposition 5.2.1 [Jeong2017]_ states that all rearrangement
        cones are proper when ``p`` is less than ``ambient_dim``::

            sage: all( cones.rearrangement(p,ambient_dim).is_proper()
            ....:              for ambient_dim in range(10)
            ....:              for p in range(1, ambient_dim) )
            True

        Jeong's Corollary 5.2.4 [Jeong2017]_ states that if ``p`` is
        either one or ``ambient_dim``, then the Lyapunov rank of the
        rearrangement cone in ``ambient_dim`` dimensions is
        ``ambient_dim``. Moreover for all other values of ``p``, its
        Lyapunov rank is one::

            sage: all( cones.rearrangement(p,ambient_dim).lyapunov_rank()
            ....:      ==
            ....:      ambient_dim
            ....:              for ambient_dim in range(2, 10)
            ....:              for p in [1, ambient_dim-1] )
            True
            sage: all( cones.rearrangement(p,ambient_dim).lyapunov_rank() == 1
            ....:              for ambient_dim in range(3, 10)
            ....:              for p in range(2, ambient_dim-1) )
            True

        TESTS:

        Jeong's Proposition 5.2.1 [Jeong2017]_ states that rearrangement
        cones are permutation-invariant::

            sage: ambient_dim = ZZ.random_element(2,10).abs()
            sage: p = ZZ.random_element(1,ambient_dim)
            sage: K = cones.rearrangement(p,ambient_dim)
            sage: P = SymmetricGroup(ambient_dim).random_element().matrix()
            sage: all( K.contains(P*r) for r in K )
            True

        The smallest ``p`` components of every element of the rearrangement
        cone should sum to a nonnegative number. In other words, the
        generators really are what we think they are::

            sage: set_random_seed()
            sage: def _has_rearrangement_property(v,p):
            ....:     return sum( sorted(v)[0:p] ) >= 0
            sage: all( _has_rearrangement_property(
            ....:      cones.rearrangement(p,ambient_dim).random_element(),
            ....:      p
            ....:    )
            ....:    for ambient_dim in range(2, 10)
            ....:    for p in range(1, ambient_dim-1)
            ....: )
            True

        Jeong's Proposition 5.2.1 [Jeong2017]_ states that the rearrangenent
        cone of order ``p`` is contained in the rearrangement cone of
        order ``p + 1``::

            sage: set_random_seed()
            sage: ambient_dim = ZZ.random_element(2,10)
            sage: p = ZZ.random_element(1,ambient_dim)
            sage: K1 = cones.rearrangement(p,ambient_dim)
            sage: K2 = cones.rearrangement(p+1,ambient_dim)
            sage: all( x in K2 for x in K1 )
            True

        Jeong's Proposition 5.2.1 [Jeong2017]_ states that the rearrangement
        cone of order ``p`` is linearly isomorphic to the rearrangement
        cone of order ``ambient_dim - p`` when ``p`` is less than
        ``ambient_dim``::

            sage: set_random_seed()
            sage: ambient_dim = ZZ.random_element(2,10)
            sage: p = ZZ.random_element(1,ambient_dim)
            sage: K1 = cones.rearrangement(p,ambient_dim)
            sage: K2 = cones.rearrangement(ambient_dim-p, ambient_dim)
            sage: Mp = ((1/p)*matrix.ones(QQ,ambient_dim)
            ....:    - matrix.identity(QQ,ambient_dim))
            sage: Cone( (Mp*K2.rays()).columns() ).is_equivalent(K1)
            True

        The order ``p`` should be between ``1`` and ``ambient_dim``,
        inclusive::

            sage: cones.rearrangement(0,3)
            Traceback (most recent call last):
            ...
            ValueError: order p=0 should be between 1 and ambient_dim=3,
            inclusive
            sage: cones.rearrangement(5,3)
            Traceback (most recent call last):
            ...
            ValueError: order p=5 should be between 1 and ambient_dim=3,
            inclusive

        If a ``lattice`` was given, it is actually used::

            sage: L = ToricLattice(3, 'M')
            sage: cones.rearrangement(2, 3, lattice=L)
            3-d cone in 3-d lattice M

        Unless the rank of the lattice disagrees with ``ambient_dim``::

            sage: L = ToricLattice(1, 'M')
            sage: cones.rearrangement(2, 3, lattice=L)
            Traceback (most recent call last):
            ...
            ValueError: lattice rank=1 and ambient_dim=3 are incompatible

        We also get an error if neither the ambient dimension nor lattice
        are specified::

            sage: cones.rearrangement(3)
            Traceback (most recent call last):
            ...
            ValueError: either the ambient dimension or the lattice must
            be specified

        """
        if ambient_dim is None and lattice is None:
            raise ValueError('either the ambient dimension or the lattice '
                             'must be specified')

        if ambient_dim is None:
            ambient_dim = lattice.rank()

        if lattice is None:
            lattice = ToricLattice(ambient_dim)

        if lattice.rank() != ambient_dim:
            raise ValueError('lattice rank=%d and ambient_dim=%d '
                             'are incompatible' % (lattice.rank(), ambient_dim))

        if p < 1 or p > ambient_dim:
            raise ValueError('order p=%d should be between 1 '
                             'and ambient_dim=%d, inclusive' % (p,ambient_dim))

        I = matrix.identity(ZZ,ambient_dim)
        M = matrix.ones(ZZ,ambient_dim) - p*I
        G = matrix.identity(ZZ,ambient_dim).rows() + M.rows()
        return Cone(G, lattice=lattice)


    def schur(self, ambient_dim=None, lattice=None):
        r"""
        The Schur cone in ``ambient_dim`` dimensions, or living
        in ``lattice``.

        INPUT:

        - ``ambient_dim`` -- (default: ``None``) the dimension of the
          ambient space, a nonnegative integer.

        - ``lattice`` -- (default: ``None``) the toric lattice in which
          the cone will live.

        If ``ambient_dim`` is omitted, then it will be inferred from the
        rank of ``lattice``. If the ``lattice`` is omitted, then the
        default lattice of rank ``ambient_dim`` will be used. A
        ``ValueError`` is raised if neither ``ambient_dim`` nor
        ``lattice`` are specified.

        It is a ``ValueError`` to specify both ``ambient_dim`` and ``lattice``
        unless the rank of ``lattice`` is equal to ``ambient_dim``.

        OUTPUT:

        A :class:`.ConvexRationalPolyhedralCone` representing the Schur
        cone living in ``lattice``, with ambient dimension ``ambient_dim``.
        Each generating ray has the integer ring as its base ring.

        A ``ValueError`` is raised if,

        - Neither the ambient dimension nor the ambient lattice are specified.

        - Both the ambient dimension and the ambient lattice are specified,
          but they are incompatible.

        REFERENCES:

        - [GS2010]_

        - [IS2005]_

        - [SS2016]_

        EXAMPLES:

        Verify the claim [SS2016]_ that the maximal angle between any two
        generators of the Schur cone and the nonnegative orthant in
        dimension five is `\left(3/4\right)\pi`::

            sage: P = cones.schur(5)
            sage: Q = cones.nonnegative_orthant(5)
            sage: G = ( g.change_ring(QQbar).normalized() for g in P )
            sage: H = ( h.change_ring(QQbar).normalized() for h in Q )
            sage: actual = max(arccos(u.inner_product(v)) for u in G for v in H)
            sage: expected = 3*pi/4
            sage: abs(actual - expected).n() < 1e-12
            True

        The dual of the Schur cone is the "downward monotonic cone"
        [GS2010]_, whose elements' entries are in non-increasing order::

            sage: set_random_seed()
            sage: ambient_dim = ZZ.random_element(10)
            sage: K = cones.schur(ambient_dim).dual()
            sage: x = K.random_element()
            sage: all( x[i] >= x[i+1] for i in range(ambient_dim-1) )
            True

        TESTS:

        We get the trivial cone when ``ambient_dim`` is zero::

            sage: cones.schur(0).is_trivial()
            True

        The Schur cone induces the majorization ordering, as in Iusem
        and Seeger's [IS2005]_ Example 7.3::

            sage: set_random_seed()
            sage: def majorized_by(x,y):
            ....:     return (all(sum(x[0:i]) <= sum(y[0:i])
            ....:                 for i in range(x.degree()-1))
            ....:             and sum(x) == sum(y))
            sage: ambient_dim = ZZ.random_element(10)
            sage: V = VectorSpace(QQ, ambient_dim)
            sage: S = cones.schur(ambient_dim)
            sage: majorized_by(V.zero(), S.random_element())
            True
            sage: x = V.random_element()
            sage: y = V.random_element()
            sage: majorized_by(x,y) == ( (y-x) in S )
            True

        If a ``lattice`` was given, it is actually used::

            sage: L = ToricLattice(3, 'M')
            sage: cones.schur(3, lattice=L)
            2-d cone in 3-d lattice M

        Unless the rank of the lattice disagrees with ``ambient_dim``::

            sage: L = ToricLattice(1, 'M')
            sage: cones.schur(3, lattice=L)
            Traceback (most recent call last):
            ...
            ValueError: lattice rank=1 and ambient_dim=3 are incompatible

        We also get an error if no arguments are given::

            sage: cones.schur()
            Traceback (most recent call last):
            ...
            ValueError: either the ambient dimension or the lattice must
            be specified

        """
        if ambient_dim is None and lattice is None:
            raise ValueError('either the ambient dimension or the lattice '
                             'must be specified')

        if ambient_dim is None:
            ambient_dim = lattice.rank()

        if lattice is None:
            lattice = ToricLattice(ambient_dim)

        if lattice.rank() != ambient_dim:
            raise ValueError('lattice rank=%d and ambient_dim=%d '
                             'are incompatible' % (lattice.rank(), ambient_dim))

        def _f(i,j):
            if i == j:
                return 1
            elif j - i == 1:
                return -1
            else:
                return 0

        # The "max" below catches the trivial case where ambient_dim == 0.
        S = matrix(ZZ, max(0,ambient_dim-1), ambient_dim, _f)

        return Cone(S.rows(), lattice)


    def trivial(self, ambient_dim=None, lattice=None):
        r"""
        The trivial cone with no nonzero generators in ``ambient_dim``
        dimensions, or living in ``lattice``.

        INPUT:

        - ``ambient_dim`` -- (default: ``None``) the dimension of the
          ambient space, a nonnegative integer.

        - ``lattice`` -- (default: ``None``) the toric lattice in which
          the cone will live.

        If ``ambient_dim`` is omitted, then it will be inferred from the
        rank of ``lattice``. If the ``lattice`` is omitted, then the
        default lattice of rank ``ambient_dim`` will be used. A
        ``ValueError`` is raised if neither ``ambient_dim`` nor
        ``lattice`` are specified.

        It is a ``ValueError`` to specify both ``ambient_dim`` and ``lattice``
        unless the rank of ``lattice`` is equal to ``ambient_dim``.

        OUTPUT:

        A :class:`.ConvexRationalPolyhedralCone` representing the
        trivial cone with no nonzero generators living in ``lattice``,
        with ambient dimension ``ambient_dim``.

        A ``ValueError`` is raised if,

        - Neither the ambient dimension nor the ambient lattice are specified.

        - Both the ambient dimension and the ambient lattice are specified,
          but they are incompatible.

        EXAMPLES:

        Construct the trivial cone, containing only the origin, in three
        dimensions::

            sage: cones.trivial(3)
            0-d cone in 3-d lattice N

        If a ``lattice`` is given, the trivial cone will live in that
        lattice::

            sage: L = ToricLattice(3, 'M')
            sage: cones.trivial(3, lattice=L)
            0-d cone in 3-d lattice M

        TESTS:

        We can construct the trivial cone in a trivial ambient space::

            sage: cones.trivial(0)
            0-d cone in 0-d lattice N

        An error is raised if the rank of the lattice disagrees with
        ``ambient_dim``::

            sage: L = ToricLattice(1, 'M')
            sage: cones.trivial(3, lattice=L)
            Traceback (most recent call last):
            ...
            ValueError: lattice rank=1 and ambient_dim=3 are incompatible

        We also get an error if no arguments are given::

            sage: cones.trivial()
            Traceback (most recent call last):
            ...
            ValueError: either the ambient dimension or the lattice must
            be specified

        """
        if ambient_dim is None and lattice is None:
            raise ValueError('either the ambient dimension or the lattice '
                             'must be specified')

        if ambient_dim is None:
            ambient_dim = lattice.rank()

        if lattice is None:
            lattice = ToricLattice(ambient_dim)

        if lattice.rank() != ambient_dim:
            raise ValueError('lattice rank=%d and ambient_dim=%d '
                             'are incompatible' % (lattice.rank(), ambient_dim))

        return Cone([], lattice)


cones = ConvexRationalPolyhedralConeFactory()
