r"""
Yangians

AUTHORS:

- Travis Scrimshaw (2013-10-08): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule
from sage.algebras.associated_graded import AssociatedGradedAlgebra

import itertools

class GeneratorIndexingSet(UniqueRepresentation):
    """
    Helper class for the indexing set of the generators.
    """
    def __init__(self, index_set, level=None):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
        """
        self._index_set = index_set
        self._level = level

    def __repr__(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: GeneratorIndexingSet((1,2))
            Cartesian product of Positive integers, (1, 2), (1, 2)
            sage: GeneratorIndexingSet((1,2), 4)
            Cartesian product of (1, 2, 3, 4), (1, 2), (1, 2)
        """
        if self._level is None:
            L = PositiveIntegers()
        else:
            L = tuple(range(1, self._level+1))
        return "Cartesian product of {L}, {I}, {I}".format(L=L, I=self._index_set)

    def an_element(self):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
            sage: I.an_element()
            (3, 1, 1)
            sage: I = GeneratorIndexingSet((1,2), 5)
            sage: I.an_element()
            (3, 1, 1)
            sage: I = GeneratorIndexingSet((1,2), 1)
            sage: I.an_element()
            (1, 1, 1)
        """
        if self._level is not None and self._level < 3:
            return (1, self._index_set[0], self._index_set[0])
        return (3, self._index_set[0], self._index_set[0])

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
            sage: I.cardinality()
            +Infinity
            sage: I = GeneratorIndexingSet((1,2), level=3)
            sage: I.cardinality() == 3 * 2 * 2
            True
        """
        if self._level is not None:
            return self._level * len(self._index_set)**2
        return infinity

    __len__ = cardinality

    def __call__(self, x):
        """
        Call ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
            sage: I([1, 2])
            (1, 2)
        """
        return tuple(x)

    def __contains__(self, x):
        """
        Check containment of ``x`` in ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
            sage: (4, 1, 2) in I
            True
            sage: [4, 2, 1] in I
            True
            sage: (-1, 1, 1) in I
            False
            sage: (1, 3, 1) in I
            False

        ::

            sage: I3 = GeneratorIndexingSet((1,2), 3)
            sage: (1, 1, 2) in I3
            True
            sage: (3, 1, 1) in I3
            True
            sage: (4, 1, 1) in I3
            False
        """
        return (isinstance(x, (tuple, list)) and len(x) == 3
                and x[0] in ZZ and x[0] > 0
                and (self._level is None or x[0] <= self._level)
                and x[1] in self._index_set
                and x[2] in self._index_set)

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: from sage.algebras.yangian import GeneratorIndexingSet
            sage: I = GeneratorIndexingSet((1,2))
            sage: it = iter(I)
            sage: [it.next() for dummy in range(5)]
            [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1)]

            sage: I = GeneratorIndexingSet((1,2), 3)
            sage: list(I)
            [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2),
             (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2),
             (3, 1, 1), (3, 1, 2), (3, 2, 1), (3, 2, 2)]
        """
        I = self._index_set
        if self._level is not None:
            for x in itertools.product(range(1, self._level+1), I, I):
                yield x
            return
        for i in PositiveIntegers():
            for x in itertools.product(I, I):
                yield (i, x[0], x[1])

class Yangian(CombinatorialFreeModule):
    r"""
    The Yangian `Y(\mathfrak{gl}_n)`.

    Let `A` be a commutative ring with unity.  The *Yangian*
    `Y(\mathfrak{gl}_n)`, associated with the Lie algebra `\mathfrak{gl}_n`
    for `n \geq 1`, is defined to be the unital associative algebra
    generated by `\{t_{ij}^{(r)} \mid 1 \leq i,j \leq n , r \geq 1\}`
    subject to the relations

    .. MATH::

        [t_{ij}^{(M+1)}, t_{k\ell}^{(L)}] - [t_{ij}^{(M)}, t_{k\ell}^{(L+1)}]
        = t_{kj}^{(M)} t_{i\ell}^{(L)} - t_{kj}^{(L)} t_{i\ell}^{(M)},

    where `L,M \geq 0` and `t_{ij}^{(0)} = \delta_{ij} \cdot 1`. This
    system of quadratic relations is equivalent to the system of
    commutation relations

    .. MATH::

        [t_{ij}^{(r)}, t_{k\ell}^{(s)}] =
        \sum_{p=0}^{\min\{r,s\}-1} \bigl(t_{kj}^{(p)} t_{i\ell}^{(r+s-1-p)}
        - t_{kj}^{(r+s-1-p)} t_{i\ell}^{(p)} \bigr),

    where `1 \leq i,j,k,\ell \leq n` and `r,s \geq 1`.

    Let `u` be a formal variable and, for
    `1 \leq i,j \leq n`, define

    .. MATH::

        t_{ij}(u) = \delta_{ij} + \sum_{r=1}^\infty t_{ij}^{(r)} u^{-r}
        \in Y(\mathfrak{gl}_n)[\![u^{-1}]\!].

    Thus, we can write the defining relations as

    .. MATH::

        \begin{aligned}
        (u - v)[t_{ij}(u), t_{k\ell}(v)] & = t_{kj}(u) t_{i\ell}(v)
        - t_{kj}(v) t_{i\ell}(u).
        \end{aligned}

    These series can be combined into a single matrix:

    .. MATH::

        T(u) := \sum_{i,j=1}^n t_{ij}(u) \otimes E_{ij} \in Y(\mathfrak{gl}_n)
        [\![u^{-1}]\!] \otimes \operatorname{End}(\CC^n),

    where `E_{ij}` is the matrix with a `1` in the `(i,j)` position
    and zeros elsewhere.

    For `m \geq 2`, define formal variables `u_1, \ldots, u_m`.
    For any `1 \leq k \leq m`, set

    .. MATH::

        T_k(u_k) := \sum_{i,j=1}^n t_{ij}(u_k) \otimes (E_{ij})_k \in
        Y(\mathfrak{gl}_n)[\![u_1^{-1},\dots,u_m^{-1}]\!] \otimes
        \operatorname{End}(\CC^n)^{\otimes m},

    where `(E_{ij})_k = 1^{\otimes (k-1)} \otimes E_{ij} \otimes
    1^{\otimes (m-k)}`. If we consider `m = 2`, we can then also write
    the defining relations as

    .. MATH::

        R(u - v) T_1(u) T_2(v) = T_2(v) T_1(u) R(u - v),

    where `R(u) = 1 - Pu^{-1}` and `P` is the permutation operator that
    swaps the two factors. Moreover, we can write the Hopf algebra
    structure as

    .. MATH::

        \Delta \colon T(u) \mapsto T_{[1]}(u) T_{[2]}(u),
        \qquad
        S \colon T(u) \mapsto T^{-1}(u),
        \qquad
        \epsilon \colon T(u) \mapsto 1,

    where `T_{[a]} = \sum_{i,j=1}^n (1^{\otimes a-1} \otimes t_{ij}(u)
    \otimes 1^{2-a}) \otimes (E_{ij})_1`.

    We can also impose two filtrations on `Y(\mathfrak{gl}_n)`: the
    *natural* filtration `\deg t_{ij}^{(r)} = r` and the *loop*
    filtration `\deg t_{ij}^{(r)} = r - 1`. The natural filtration has
    a graded homomorphism with `U(\mathfrak{gl}_n)` by
    `t_{ij}^{(r)} \mapsto (E^r)_{ij}` and an associated graded algebra
    being polynomial algebra. Moreover, this shows a PBW theorem for
    the Yangian, that for any fixed order, we can write elements as
    unique linear combinations of ordered monomials using `t_{ij}^{(r)}`.
    For the loop filtration, the associated graded algebra is isomorphic
    (as Hopf algebras) to `U(\mathfrak{gl}_n[z])` given by
    `\overline{t}_{ij}^{(r)} \mapsto E_{ij} x^{r-1}`, where
    `\overline{t}_{ij}^{(r)}` is the image of `t_{ij}^{(r)}` in the
    `(r - 1)`-th component of `\operatorname{gr}Y(\mathfrak{gl}_n)`.

    INPUT:

    - ``base_ring`` -- the base ring
    - ``n`` -- the size `n`
    - ``level`` -- (optional) the level of the Yangian
    - ``variable_name`` -- (default: ``'t'``) the name of the variable
    - ``filtration`` -- (default: ``'loop'``) the filtration and can be
      one of the following:

      * ``'natural'`` -- the filtration is given by `\deg t_{ij}^{(r)} = r`
      * ``'loop'`` -- the filtration is given by `\deg t_{ij}^{(r)} = r - 1`

    .. TODO::

        Implement the antipode.

    EXAMPLES::

        sage: Y = Yangian(QQ, 4)
        sage: t = Y.algebra_generators()
        sage: t[6,2,1] * t[2,3,2]
        -t(1)[2,2]*t(6)[3,1] + t(1)[3,1]*t(6)[2,2]
         + t(2)[3,2]*t(6)[2,1] - t(7)[3,1]
        sage: t[6,2,1] * t[3,1,4]
        t(1)[1,1]*t(7)[2,4] + t(1)[1,4]*t(6)[2,1] - t(1)[2,1]*t(6)[1,4]
         - t(1)[2,4]*t(7)[1,1] + t(2)[1,1]*t(6)[2,4] - t(2)[2,4]*t(6)[1,1]
         + t(3)[1,4]*t(6)[2,1] + t(6)[2,4] + t(8)[2,4]

    We check that the natural filtration has a homomorphism
    to `U(\mathfrak{gl}_n)` as algebras::

        sage: Y = Yangian(QQ, 4, filtration='natural')
        sage: t = Y.algebra_generators()
        sage: gl4 = lie_algebras.gl(QQ, 4)
        sage: Ugl4 = gl4.pbw_basis()
        sage: E = matrix(Ugl4, 4, 4, Ugl4.gens())
        sage: Esq = E^2
        sage: t[2,1,3] * t[1,2,1]
        t(1)[2,1]*t(2)[1,3] - t(2)[2,3]
        sage: Esq[0,2] * E[1,0] == E[1,0] * Esq[0,2] - Esq[1,2]
        True

        sage: Em = [E^k for k in range(1,5)]
        sage: S = list(t.some_elements())[:30:3]
        sage: def convert(x):
        ....:     return sum(c * prod(Em[t[0]-1][t[1]-1,t[2]-1] ** e
        ....:                         for t,e in m._sorted_items())
        ....:                for m,c in x)
        sage: for x in S:
        ....:     for y in S:
        ....:         ret = x * y
        ....:         rhs = convert(x) * convert(y)
        ....:         assert rhs == convert(ret)
        ....:         assert ret.maximal_degree() == rhs.maximal_degree()

    REFERENCES:

    - :wikipedia:`Yangian`
    - [MNO1994]_
    - [Mol2007]_
    """
    @staticmethod
    def __classcall_private__(cls, base_ring, n, level=None,
                              variable_name='t', filtration='loop'):
        """
        Return the correct parent based upon input.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y2 = Yangian(QQ, 4)
            sage: Y is Y2
            True
            sage: YL = Yangian(QQ, 4, 3)
            sage: YL2 = Yangian(QQ, 4, 3)
            sage: YL is YL2
            True
        """
        if filtration not in ['natural', 'loop']:
            raise ValueError("invalid filtration")

        if level is not None:
            return YangianLevel(base_ring, n, level, variable_name, filtration)
        # We need to specify the parameter name for pickling, so it doesn't pass
        #   ``variable_name`` as ``level``
        return super(Yangian, cls).__classcall__(cls, base_ring, n,
                                                 variable_name=variable_name,
                                                 filtration=filtration)

    def __init__(self, base_ring, n, variable_name, filtration):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: TestSuite(Y).run(skip="_test_antipode") # Not implemented
        """
        self._n = n
        self._filtration = filtration
        category = HopfAlgebrasWithBasis(base_ring).Filtered()
        self._index_set = tuple(range(1,n+1))
        # The keys for the basis are tuples (l, i, j)
        indices = GeneratorIndexingSet(self._index_set)
        # We note that the generators are non-commutative, but we always sort
        #   them, so they are, in effect, indexed by the free abelian monoid
        basis_keys = IndexedFreeAbelianMonoid(indices, bracket=False,
                                              prefix=variable_name)
        CombinatorialFreeModule.__init__(self, base_ring, basis_keys,
                                         sorting_key=Yangian._term_key,
                                         prefix=variable_name, category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Yangian(QQ, 4)
            Yangian of gl(4) in the loop filtration over Rational Field
            sage: Yangian(QQ, 4, filtration='natural')
            Yangian of gl(4) in the natural filtration over Rational Field
        """
        return "Yangian of gl({}) in the {} filtration over {}".format(self._n, self._filtration, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(Yangian(QQ, 4))
            Y(\mathfrak{gl}_{4}, \Bold{Q})
        """
        from sage.misc.latex import latex
        return "Y(\\mathfrak{{gl}}_{{{}}}, {})".format(self._n, latex(self.base_ring()))

    @staticmethod
    def _term_key(x):
        """
        Compute a key for ``x`` for comparisons.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: x = Y.gen(2, 1, 1).leading_support()
            sage: Yangian._term_key(x)
            (-1, [((2, 1, 1), 1)])
        """
        return (-len(x), x._sorted_items())

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: I = Y._indices
            sage: Y._repr_term(I.gen((3,1,2))^2 * I.gen((4,3,1)))
            't(3)[1,2]^2*t(4)[3,1]'
            sage: Y._repr_term(Y.one_basis())
            '1'
        """
        if len(m) == 0:
            return '1'
        prefix = self.prefix()
        return '*'.join(prefix + '({})[{},{}]'.format(r,i,j)
                        + ('^{}'.format(exp) if exp > 1 else '')
                        for (r,i,j), exp in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: I = Y._indices
            sage: Y._latex_term(I.gen((3,1,2))^2 * I.gen((4,3,1)))
            '\\left(t^{(3)}_{1,2}\\right)^{2} t^{(4)}_{3,1}'
            sage: Y._latex_term(Y.one_basis())
            '1'
        """
        if len(m) == 0:
            return '1'

        prefix = self.prefix()
        def term(r, i, j, exp):
            s = prefix + '^{{({})}}_{{{},{}}}'.format(r,i,j)
            if exp == 1:
                return s
            return '\\left({}\\right)^{{{}}}'.format(s, exp)
        return ' '.join(term(r, i, j, exp) for (r,i,j), exp in m._sorted_items())

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Yn = Yangian(QQ, 4, filtration='natural')
            sage: Yn(Y.an_element()) == Yn.an_element()
            True
            sage: Y(Yn.an_element()) == Y.an_element()
            True
            sage: Y6 = Yangian(QQ, 4, level=6, filtration='natural')
            sage: Y(Y6.an_element())
            t(1)[1,1]*t(1)[1,2]^2*t(1)[1,3]^3*t(3)[1,1]
        """
        if isinstance(x, CombinatorialFreeModule.Element):
            if isinstance(x.parent(), Yangian) and x.parent()._n <= self._n:
                R = self.base_ring()
                return self._from_dict({i: R(c) for i,c in x}, coerce=False)
        return super(Yangian, self)._element_constructor_(x)

    def gen(self, r, i=None, j=None):
        """
        Return the generator `t^{(r)}_{ij}` of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2, 1, 3)
            t(2)[1,3]
            sage: Y.gen(12, 2, 1)
            t(12)[2,1]
            sage: Y.gen(0, 1, 1)
            1
            sage: Y.gen(0, 1, 3)
            0
        """
        if i is None and j is None:
            r,i,j = r
        if r == 0:
            if i == j:
                return self.one()
            return self.zero()
        m = self._indices.gen((r,i,j))
        return self.element_class(self, {m: self.base_ring().one()})

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.algebra_generators()
            Lazy family (generator(i))_{i in Cartesian product of
             Positive integers, (1, 2, 3, 4), (1, 2, 3, 4)}
        """
        return Family(self._indices._indices, self.gen, name="generator")

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.one_basis()
            1
        """
        return self._indices.one()

    def degree_on_basis(self, m):
        """
        Return the degree of the monomial index by ``m``.

        The degree of `t_{ij}^{(r)}` is equal to `r - 1` if ``filtration =
        'loop'`` and is equal to `r` if ``filtration = 'natural'``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.degree_on_basis(Y.gen(2,1,1).leading_support())
            1
            sage: x = Y.gen(5,2,3)^4
            sage: Y.degree_on_basis(x.leading_support())
            16
            sage: elt = Y.gen(10,3,1) * Y.gen(2,1,1) * Y.gen(1,2,4); elt
            t(1)[1,1]*t(1)[2,4]*t(10)[3,1] - t(1)[2,4]*t(1)[3,1]*t(10)[1,1]
             + t(1)[2,4]*t(2)[1,1]*t(10)[3,1] + t(1)[2,4]*t(10)[3,1]
             + t(1)[2,4]*t(11)[3,1]
            sage: for s in sorted(elt.support(), key=str): s, Y.degree_on_basis(s)
            (t(1, 1, 1)*t(1, 2, 4)*t(10, 3, 1), 9)
            (t(1, 2, 4)*t(1, 3, 1)*t(10, 1, 1), 9)
            (t(1, 2, 4)*t(10, 3, 1), 9)
            (t(1, 2, 4)*t(11, 3, 1), 10)
            (t(1, 2, 4)*t(2, 1, 1)*t(10, 3, 1), 10)

            sage: Y = Yangian(QQ, 4, filtration='natural')
            sage: Y.degree_on_basis(Y.gen(2,1,1).leading_support())
            2
            sage: x = Y.gen(5,2,3)^4
            sage: Y.degree_on_basis(x.leading_support())
            20
            sage: elt = Y.gen(10,3,1) * Y.gen(2,1,1) * Y.gen(1,2,4)
            sage: for s in sorted(elt.support(), key=str): s, Y.degree_on_basis(s)
            (t(1, 1, 1)*t(1, 2, 4)*t(10, 3, 1), 12)
            (t(1, 2, 4)*t(1, 3, 1)*t(10, 1, 1), 12)
            (t(1, 2, 4)*t(10, 3, 1), 11)
            (t(1, 2, 4)*t(11, 3, 1), 12)
            (t(1, 2, 4)*t(2, 1, 1)*t(10, 3, 1), 13)
        """
        if self._filtration == 'natural':
            return sum(r[0][0] * r[1] for r in m._monomial.items())
        return sum(max(0, r[0][0] - 1) * r[1] for r in m._monomial.items())

    def graded_algebra(self):
        """
        Return the associated graded algebra of ``self``.

        EXAMPLES::

            sage: Yangian(QQ, 4).graded_algebra()
            Graded Algebra of Yangian of gl(4) in the loop filtration over Rational Field
            sage: Yangian(QQ, 4, filtration='natural').graded_algebra()
            Graded Algebra of Yangian of gl(4) in the natural filtration over Rational Field
        """
        if self._filtration == 'natural':
            return GradedYangianNatural(self)
        return GradedYangianLoop(self)

    def dimension(self):
        """
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.dimension()
            +Infinity
        """
        return infinity

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of two monomials given by ``x`` and ``y``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(12, 2, 1) * Y.gen(2, 1, 1) # indirect doctest
            t(1)[1,1]*t(12)[2,1] - t(1)[2,1]*t(12)[1,1]
             + t(2)[1,1]*t(12)[2,1] + t(12)[2,1] + t(13)[2,1]
        """
        # If x or y indexed by the identity element, it is 1, so return the other
        if len(x) == 0:
            return self.monomial(y)
        if len(y) == 0:
            return self.monomial(x)
        # If it's smaller, just add it to the front
        if x.trailing_support() <= y.leading_support():
            return self.monomial(x * y)

        # The computation is done on generators, so apply generators one at
        #   a time until all have been applied
        if len(x) != 1:
            I = self._indices
            cur = self.monomial(y)
            for gen,exp in reversed(x._sorted_items()):
                for i in range(exp):
                    cur = self.monomial(I.gen(gen)) * cur
            return cur

        # If we are both generators, then apply the basic computation
        if len(y) == 1:
            return self.product_on_gens(tuple(x.support()[0]), tuple(y.support()[0]))

        # Otherwise we need to commute it along
        I = self._indices
        cur = self.monomial(x)
        for gen,exp in y._sorted_items():
            for i in range(exp):
                cur = cur * self.monomial(I.gen(gen))
        return cur

    @cached_method
    def product_on_gens(self, a, b):
        r"""
        Return the product on two generators indexed by ``a`` and ``b``.

        We assume `(r, i, j) \geq (s, k, \ell)`, and we start with the basic
        relation:

        .. MATH::

            [t_{ij}^{(r)}, t_{k\ell}^{(s)}] - [t_{ij}^{(r-1)}, t_{k\ell}^{(s+1)}]
            = t_{kj}^{(r-1)} t_{i\ell}^{(s)} - t_{kj}^{(s)} t_{i\ell}^{(r-1)}.

        Solving for the first term and using induction we get:

        .. MATH::

            [t_{ij}^{(r)}, t_{k\ell}^{(s)}] = \sum_{a=1}^s \left(
            t_{kj}^{(a-1)} t_{i\ell}^{(r+s-a)} - t_{kj}^{(r+s-a)}
            t_{i\ell}^{(a-1)} \right).

        Next applying induction on this we get

        .. MATH::

            t_{ij}^{(r)} t_{k\ell}^{(s)} = t_{k\ell}^{(s)} t_{ij}^{(r)} +
            \sum C_{abcd}^{m\ell} t_{ab}^{(m)} t_{cd}^{(\ell)}

        where `m + \ell < r + s` and `t_{ab}^{(m)} < t_{cd}^{(\ell)}`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.product_on_gens((2,1,1), (12,2,1))
            t(2)[1,1]*t(12)[2,1]
            sage: Y.gen(2, 1, 1) * Y.gen(12, 2, 1)
            t(2)[1,1]*t(12)[2,1]
            sage: Y.product_on_gens((12,2,1), (2,1,1))
            t(1)[1,1]*t(12)[2,1] - t(1)[2,1]*t(12)[1,1]
             + t(2)[1,1]*t(12)[2,1] + t(12)[2,1] + t(13)[2,1]
            sage: Y.gen(12, 2, 1) * Y.gen(2, 1, 1)
            t(1)[1,1]*t(12)[2,1] - t(1)[2,1]*t(12)[1,1]
             + t(2)[1,1]*t(12)[2,1] + t(12)[2,1] + t(13)[2,1]
        """
        I = self._indices
        if a <= b:
            return self.monomial(I.gen(a) * I.gen(b))

        # This is the special term of x = 1
        x1 = self.zero()
        if b[1] == a[2]:
            x1 += self.monomial( I.gen((a[0]+b[0]-1, a[1], b[2])) )
        if a[1] == b[2]:
            x1 -= self.monomial( I.gen((a[0]+b[0]-1, b[1], a[2])) )

        return self.monomial(I.gen(b) * I.gen(a)) + x1 + self.sum(
                self.monomial( I.gen((x-1, b[1], a[2])) * I.gen((a[0]+b[0]-x, a[1], b[2])) )
                - self.product_on_gens( (a[0]+b[0]-x, b[1], a[2]), (x-1, a[1], b[2]) )
                for x in range(2, b[0]+1))

    def coproduct_on_basis(self, m):
        """
        Return the coproduct on the basis element indexed by ``m``.

        The coproduct `\Delta\colon Y(\mathfrak{gl}_n) \longrightarrow
        Y(\mathfrak{gl}_n) \otimes Y(\mathfrak{gl}_n)` is defined by

        .. MATH::

            \Delta(t_{ij}(u)) = \sum_{a=1}^n t_{ia}(u) \otimes t_{aj}(u).

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2,1,1).coproduct() # indirect doctest
            1 # t(2)[1,1] + t(1)[1,1] # t(1)[1,1] + t(1)[1,2] # t(1)[2,1]
             + t(1)[1,3] # t(1)[3,1] + t(1)[1,4] # t(1)[4,1] + t(2)[1,1] # 1
            sage: Y.gen(2,3,1).coproduct()
            1 # t(2)[3,1] + t(1)[3,1] # t(1)[1,1] + t(1)[3,2] # t(1)[2,1]
             + t(1)[3,3] # t(1)[3,1] + t(1)[3,4] # t(1)[4,1] + t(2)[3,1] # 1
            sage: Y.gen(2,2,3).coproduct()
            1 # t(2)[2,3] + t(1)[2,1] # t(1)[1,3] + t(1)[2,2] # t(1)[2,3]
             + t(1)[2,3] # t(1)[3,3] + t(1)[2,4] # t(1)[4,3] + t(2)[2,3] # 1
        """
        T = self.tensor_square()
        I = self._indices
        return T.prod(T.monomial( (I.one(), I.gen((a[0],a[1],a[2]))) )
                      + T.monomial( (I.gen((a[0],a[1],a[2])), I.one()) )
                      + T.sum_of_terms([(( I.gen((s,a[1],k)), I.gen((a[0]-s,k,a[2])) ), 1)
                                        for k in range(1, self._n+1)
                                        for s in range(1, a[0])])
                      for a,exp in m._sorted_items() for p in range(exp))

    def counit_on_basis(self, m):
        """
        Return the counit on the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2,3,1).counit() # indirect doctest
            0
            sage: Y.gen(0,0,0).counit()
            1
        """
        if len(m) == 0:
            return self.base_ring().one()
        return self.base_ring().zero()

class YangianLevel(Yangian):
    r"""
    The Yangian `Y_{\ell}(\mathfrak{gl_n})` of level `\ell`.

    The Yangian of level `\ell` is the quotient of the Yangian
    `Y(\mathfrak{gl}_n)` by the two-sided ideal generated by `t_{ij}^{(r)}`
    for all `r > p` and all `i,j \in \{1, \ldots, n\}`.

    EXAMPLES::

        sage: Y = Yangian(QQ, 4, 3)
        sage: elt = Y.gen(3,2,1) * Y.gen(1,1,3)
        sage: elt * Y.gen(1, 1, 2)
        t(1)[1,2]*t(1)[1,3]*t(3)[2,1] + t(1)[1,2]*t(3)[2,3]
         - t(1)[1,3]*t(3)[1,1] + t(1)[1,3]*t(3)[2,2] - t(3)[1,3]
    """
    def __init__(self, base_ring, n, level, variable_name, filtration):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: TestSuite(Y).run(skip="_test_antipode")
        """
        self._level = level
        self._n = n
        self._filtration = filtration
        category = HopfAlgebrasWithBasis(base_ring)#.Filtered() # TODO - once implemented
        self._index_set = tuple(range(1,n+1))
        # The keys for the basis are tuples (l, i, j)
        indices = GeneratorIndexingSet(self._index_set, level)
        # We note that the generators are non-commutative, but we always sort
        #   them, so they are, in effect, indexed by the free abelian monoid
        basis_keys = IndexedFreeAbelianMonoid(indices, bracket=False, prefix=variable_name)
        CombinatorialFreeModule.__init__(self, base_ring, basis_keys,
                                         prefix=variable_name, category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Yangian(QQ, 4, 3)
            Yangian of level 3 of gl(4) in the loop filtration over Rational Field
        """
        return "Yangian of level {} of gl({}) in the {} filtration over {}".format(
                        self._level, self._n, self._filtration, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(Yangian(QQ, 4, level=5))
            Y_{5}(\mathfrak{gl}_{4}, \Bold{Q})
        """
        from sage.misc.latex import latex
        return "Y_{{{}}}(\\mathfrak{{gl}}_{{{}}}, {})".format(
                        self._level, self._n, latex(self.base_ring()))

    def _coerce_map_from_(self, R):
        """
        Return ``True`` or the coercion if there exists a coerce
        map from ``R``.

        EXAMPLES::

            sage: Y5 = Yangian(QQ, 7, level=5)
            sage: Y = Yangian(QQ, 3)
            sage: Y5._coerce_map_from_(Y)
            Generic morphism:
              From: Yangian of gl(3) in the loop filtration over Rational Field
              To:   Yangian of level 5 of gl(7) in the loop filtration over Rational Field
            sage: phi = Y5.coerce_map_from(Y)
            sage: x = Y.gen(5,2,1) * Y.gen(4,3,2)
            sage: phi(x)
            -t(1)[2,2]*t(5)[3,1] + t(1)[3,1]*t(5)[2,2]
             - t(2)[2,1]*t(5)[3,2] + t(2)[3,2]*t(5)[2,1]
             - t(3)[2,2]*t(5)[3,1] + t(3)[3,1]*t(5)[2,2]
             + t(4)[3,2]*t(5)[2,1]

            sage: Y = Yangian(QQ, 10)
            sage: Y5.has_coerce_map_from(Y)
            False

            sage: Y10 = Yangian(QQ, 4, level=10)
            sage: phi = Y5.coerce_map_from(Y10); phi
            Generic morphism:
              From: Yangian of level 10 of gl(4) in the loop filtration over Rational Field
              To:   Yangian of level 5 of gl(7) in the loop filtration over Rational Field
            sage: x = Y10.gen(5,2,1) * Y10.gen(4,3,2)
            sage: phi(x)
            -t(1)[2,2]*t(5)[3,1] + t(1)[3,1]*t(5)[2,2]
             - t(2)[2,1]*t(5)[3,2] + t(2)[3,2]*t(5)[2,1]
             - t(3)[2,2]*t(5)[3,1] + t(3)[3,1]*t(5)[2,2]
             + t(4)[3,2]*t(5)[2,1]

            sage: Y = Yangian(QQ, 3, filtration='natural')
            sage: Y5.has_coerce_map_from(Y)
            False
        """
        if isinstance(R, Yangian) and R._n <= self._n and R._filtration == self._filtration:
            if isinstance(R, YangianLevel) and self._level > R._level:
                return False
            on_gens = lambda m: self.prod(self.gen(*a)**exp for a,exp in m._sorted_items())
            return R.module_morphism(on_gens, codomain=self)
        return super(YangianLevel, self)._coerce_map_from_(R)

    def level(self):
        """
        Return the level of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 3, 5)
            sage: Y.level()
            5
        """
        return self._level

    def defining_polynomial(self, i, j, u=None):
        """
        Return the defining polynomial of ``i`` and ``j``.

        The defining polynomial is given by:

        .. MATH::

            T_{ij}(u) = \delta_{ij} u^{\ell} + \sum_{k=1}^{\ell} t_{ij}^{(k)}
            u^{\ell-k}.

        EXAMPLES::

            sage: Y = Yangian(QQ, 3, 5)
            sage: Y.defining_polynomial(3, 2)
            t(1)[3,2]*u^4 + t(2)[3,2]*u^3 + t(3)[3,2]*u^2 + t(4)[3,2]*u + t(5)[3,2]
            sage: Y.defining_polynomial(1, 1)
            u^5 + t(1)[1,1]*u^4 + t(2)[1,1]*u^3 + t(3)[1,1]*u^2 + t(4)[1,1]*u + t(5)[1,1]
        """
        if u is None:
            u = PolynomialRing(self.base_ring(), 'u').gen(0)
        ell = self._level
        return sum(self.gen(k, i, j) * u**(ell-k) for k in range(ell+1))

    def quantum_determinant(self, u=None):
        """
        Return the quantum determinant of ``self``.

        The quantum determinant is defined by:

        .. MATH::

            \operatorname{qdet}(u) = \sum_{\sigma \in S_n} (-1)^{\sigma}
            \prod_{k=1}^n T_{\sigma(k),k}(u - k + 1).

        EXAMPLES::

            sage: Y = Yangian(QQ, 2, 2)
            sage: Y.quantum_determinant()
            u^4 + (-2 + t(1)[1,1] + t(1)[2,2])*u^3
             + (1 - t(1)[1,1] + t(1)[1,1]*t(1)[2,2] - t(1)[1,2]*t(1)[2,1]
                - 2*t(1)[2,2] + t(2)[1,1] + t(2)[2,2])*u^2
             + (-t(1)[1,1]*t(1)[2,2] + t(1)[1,1]*t(2)[2,2]
                + t(1)[1,2]*t(1)[2,1] - t(1)[1,2]*t(2)[2,1]
                - t(1)[2,1]*t(2)[1,2] + t(1)[2,2] + t(1)[2,2]*t(2)[1,1]
                - t(2)[1,1] - t(2)[2,2])*u
             - t(1)[1,1]*t(2)[2,2] + t(1)[1,2]*t(2)[2,1] + t(2)[1,1]*t(2)[2,2]
                - t(2)[1,2]*t(2)[2,1] + t(2)[2,2]
        """
        if u is None:
            u = PolynomialRing(self.base_ring(), 'u').gen(0)
        from sage.combinat.permutation import Permutations
        n = self._n
        return sum(p.sign() * prod(self.defining_polynomial(p[k], k+1, u - k)
                                   for k in range(n))
                   for p in Permutations(n))

    def gen(self, r, i=None, j=None):
        """
        Return the generator `t^{(r)}_{ij}` of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: Y.gen(2, 1, 3)
            t(2)[1,3]
            sage: Y.gen(12, 2, 1)
            0
            sage: Y.gen(0, 1, 1)
            1
            sage: Y.gen(0, 1, 3)
            0
        """
        if i is None and j is None:
            r,i,j = r
        if r > self._level:
            return self.zero()
        return Yangian.gen(self, r, i, j)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 2, 2)
            sage: Y.gens()
            (t(1)[1,1], t(2)[1,1], t(1)[1,2], t(2)[1,2], t(1)[2,1],
             t(2)[2,1], t(1)[2,2], t(2)[2,2])
        """
        return tuple(self.gen(r, i, j)
                     for i in range(1, self._n+1)
                     for j in range(1, self._n+1)
                     for r in range(1, self._level+1))

    @cached_method
    def product_on_gens(self, a, b):
        r"""
        Return the product on two generators indexed by ``a`` and ``b``.

        .. SEEALSO::

            :meth:`Yangian.product_on_gens()`

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: Y.gen(1,2,2) * Y.gen(2,1,3) # indirect doctest
            t(1)[2,2]*t(2)[1,3]
            sage: Y.gen(1,2,1) * Y.gen(2,1,3) # indirect doctest
            t(1)[2,1]*t(2)[1,3]
            sage: Y.gen(3,2,1) * Y.gen(1,1,3) # indirect doctest
            t(1)[1,3]*t(3)[2,1] + t(3)[2,3]
        """
        I = self._indices
        if a <= b:
            return self.monomial(I.gen(a) * I.gen(b))

        # This is the special term of x = 1
        x1 = self.zero()
        if a[0]+b[0]-1 <= self._level:
            if b[1] == a[2]:
                x1 += self.monomial( I.gen((a[0]+b[0]-1, a[1], b[2])) )
            if a[1] == b[2]:
                x1 -= self.monomial( I.gen((a[0]+b[0]-1, b[1], a[2])) )

        return self.monomial(I.gen(b) * I.gen(a)) + x1 + self.sum(
                self.monomial( I.gen((x-1, b[1], a[2])) * I.gen((a[0]+b[0]-x, a[1], b[2])) )
                - self.product_on_gens((a[0]+b[0]-x, b[1], a[2]), (x-1, a[1], b[2]))
                for x in range(2, b[0]+1) if a[0]+b[0]-x <= self._level)

#####################################################################
## Graded algebras

class GradedYangianBase(AssociatedGradedAlgebra):
    """
    Base class for graded algebras associated to a Yangian.
    """
    def _repr_term(self, m):
        """
        Return a string representation of the monomial indexed by ``m``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: I = grY._indices
            sage: grY._repr_term(I.gen((3,1,2))^2 * I.gen((4,3,1)))
            'tbar(3)[1,2]^2*tbar(4)[3,1]'
        """
        if len(m) == 0:
            return '1'
        prefix = self.prefix()
        return '*'.join(prefix + '({})[{},{}]'.format(r,i,j)
                        + ('^{}'.format(exp) if exp > 1 else '')
                        for (r,i,j), exp in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a latex representation of the monomial indexed by ``m``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: I = grY._indices
            sage: grY._latex_term(I.gen((3,1,2))^2 * I.gen((4,3,1)))
            '\\left(\\overline{t}^{(3)}_{1,2}\\right)^{2} \\overline{t}^{(4)}_{3,1}'
        """
        if len(m) == 0:
            return '1'

        prefix = "\\overline{{{}}}".format(self._A.prefix())
        def term(r, i, j, exp):
            s = prefix + '^{{({})}}_{{{},{}}}'.format(r,i,j)
            if exp == 1:
                return s
            return '\\left({}\\right)^{{{}}}'.format(s, exp)
        return ' '.join(term(r, i, j, exp) for (r,i,j), exp in m._sorted_items())

class GradedYangianNatural(GradedYangianBase):
    r"""
    The associated graded algebra corresponding to a Yangian
    `\operatorname{gr} Y(\mathfrak{gl}_n)` with the natural filtration
    of `\deg t_{ij}^{(r)} = r`.

    INPUT:

    - ``Y`` -- a Yangian with the natural filtration
    """
    def __init__(self, Y):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4, filtration='natural').graded_algebra()
            sage: TestSuite(grY).run(skip='_test_antipode')
        """
        if Y._filtration != 'natural':
            raise ValueError("the Yangian must have the natural filtration")
        cat = GradedHopfAlgebrasWithBasis(Y.base_ring()).Commutative()
        GradedYangianBase.__init__(self, Y, cat)

    def product_on_basis(self, x, y):
        """
        Return the product on basis elements given by the
        indices ``x`` and ``y``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4, filtration='natural').graded_algebra()
            sage: x = grY.gen(12, 2, 1) * grY.gen(2, 1, 1) # indirect doctest
            sage: x
            tbar(2)[1,1]*tbar(12)[2,1]
            sage: x == grY.gen(2, 1, 1) * grY.gen(12, 2, 1)
            True
        """
        return self.monomial(x*y)

class GradedYangianLoop(GradedYangianBase):
    r"""
    The associated graded algebra corresponding to a Yangian
    `\operatorname{gr} Y(\mathfrak{gl}_n)` with the filtration
    of `\deg t_{ij}^{(r)} = r - 1`.

    Using this filtration for the Yangian, the associated graded algebra
    is isomorphic to `U(\mathfrak{gl}_n[z])`, the universal enveloping
    algebra of the loop algebra of `\mathfrak{gl}_n`.

    INPUT:

    - ``Y`` -- a Yangian with the loop filtration
    """
    def __init__(self, Y):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: TestSuite(grY).run()
        """
        if Y._filtration != 'loop':
            raise ValueError("the Yangian must have the loop filtration")
        cat = GradedHopfAlgebrasWithBasis(Y.base_ring())
        GradedYangianBase.__init__(self, Y, cat)

    def antipode_on_basis(self, m):
        """
        Return the antipode on a basis element indexed by ``m``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: grY.antipode_on_basis(grY.gen(2,1,1).leading_support())
            -tbar(2)[1,1]

            sage: x = grY.an_element(); x
            tbar(1)[1,1]*tbar(1)[1,2]^2*tbar(1)[1,3]^3*tbar(3)[1,1]
            sage: grY.antipode_on_basis(x.leading_support())
            -tbar(1)[1,1]*tbar(1)[1,2]^2*tbar(1)[1,3]^3*tbar(3)[1,1]
             - 2*tbar(1)[1,1]*tbar(1)[1,2]*tbar(1)[1,3]^3*tbar(3)[1,2]
             - 3*tbar(1)[1,1]*tbar(1)[1,2]^2*tbar(1)[1,3]^2*tbar(3)[1,3]
             + 5*tbar(1)[1,2]^2*tbar(1)[1,3]^3*tbar(3)[1,1]
             + 10*tbar(1)[1,2]*tbar(1)[1,3]^3*tbar(3)[1,2]
             + 15*tbar(1)[1,2]^2*tbar(1)[1,3]^2*tbar(3)[1,3]
        """
        return self.prod( (-1)**exp * self.monomial(a**exp)
                          for a,exp in reversed(list(m)) )

    def coproduct_on_basis(self, m):
        """
        Return the coproduct on the basis element indexed by ``m``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: grY.coproduct_on_basis(grY.gen(2,1,1).leading_support())
            1 # tbar(2)[1,1] + tbar(2)[1,1] # 1
            sage: grY.gen(2,3,1).coproduct()
            1 # tbar(2)[3,1] + tbar(2)[3,1] # 1
        """
        T = self.tensor_square()
        I = self._indices
        one = I.one()
        return T.prod(T.sum_of_monomials([(one, a), (a, one)])
                      for a,exp in m for p in range(exp))

    def counit_on_basis(self, m):
        """
        Return the antipode on the basis element indexed by ``m``.

        EXAMPLES::

            sage: grY = Yangian(QQ, 4).graded_algebra()
            sage: grY.counit_on_basis(grY.gen(2,3,1).leading_support())
            0
            sage: grY.gen(0,0,0).counit()
            1
        """
        if len(m) == 0:
            return self.base_ring().one()
        return self.base_ring().zero()

