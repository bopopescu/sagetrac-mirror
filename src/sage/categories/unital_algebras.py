r"""
Unital algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.magmatic_algebras import MagmaticAlgebras

class UnitalAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of non-associative algebras over a given base ring.

    A non-associative algebra over a ring `R` is a module over `R`
    which s also a unital magma.

    .. WARNING::

        Until :trac:`15043` is implemented, :class:`Algebras` is the
        category of associative unital algebras; thus, unlike the name
        suggests, :class:`UnitalAlgebras` is not a subcategory of
        :class:`Algebras` but of
        :class:`~.magmatic_algebras.MagmaticAlgebras`.

    EXAMPLES::

        sage: from sage.categories.unital_algebras import UnitalAlgebras
        sage: C = UnitalAlgebras(ZZ); C
        Category of unital algebras over Integer Ring

    TESTS::

        sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
        sage: C is MagmaticAlgebras(ZZ).Unital()
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (MagmaticAlgebras, "Unital")

    class ParentMethods:
        def from_base_ring(self, r):
            """
            Return the canonical embedding of ``r`` into ``self``.

            INPUT:

            - ``r`` -- an element of ``self.base_ring()``

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: A.from_base_ring(1)
                B[word: ]
            """
            return self.one()._lmul_(r)

        def _coerce_map_from_(self, X):
            """
            Return a coercion map from `X` to ``self``, or ``None``.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: coercion_model = sage.structure.element.get_coercion_model()
                sage: coercion_model.discover_coercion(QQ, A)
                ((map internal to coercion system -- copy before use)
                 Generic morphism:
                  From: Rational Field
                  To:   An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field, None)
                sage: A(1)          # indirect doctest
                B[word: ]

            Check that :trac:`19225` is solved::

                    sage: A = cartesian_product((QQ['z'],)); A
                    The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)
                    sage: A.base_ring()
                    Rational Field
                    sage: A(1)
                    (1,)
            """
            base_ring = self.base_ring()
            if X is base_ring:
                from sage.categories.sets_cat import Sets
                H = Hom(base_ring, self, Sets())

                # Idea: There is a generic method "from_base_ring",
                # that just multiplies by the multiplicative unit.
                # However, the unit is constructed repeatedly, which
                # may be slow, so we store the unit.
                #
                # However, if there is a specialised from_base_ring
                # method, then it should be used!
                try:
                    has_custom_conversion = self.category().parent_class.from_base_ring.__func__ is not self.from_base_ring.__func__
                except AttributeError:
                    # Sometimes from_base_ring is a lazy attribute
                    has_custom_conversion = True
                if has_custom_conversion:
                    return SetMorphism(function=self.from_base_ring, parent=H)
                try:
                    one = self.one()
                    return SetMorphism(function=one._lmul_, parent=H)
                except (NotImplementedError, AttributeError, TypeError):
                    # It is possible that an_element or lmul are not
                    # implemented.
                    pass
            return self._coerce_map_via([base_ring], X)


    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            @abstract_method(optional = True)
            def one_basis(self):
                """
                When the one of an algebra with basis is an element of
                this basis, this optional method can return the index of
                this element. This is used to provide a default
                implementation of :meth:`.one`, and an optimized default
                implementation of :meth:`.from_base_ring`.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word: 
                    sage: A.one()
                    B[word: ]
                    sage: A.from_base_ring(4)
                    4*B[word: ]
                """

            @cached_method
            def one_from_one_basis(self):
                """
                Return the one of the algebra, as per
                :meth:`Monoids.ParentMethods.one()
                <sage.categories.monoids.Monoids.ParentMethods.one>`

                By default, this is implemented from
                :meth:`.one_basis`, if available.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word: 
                    sage: A.one_from_one_basis()
                    B[word: ]
                    sage: A.one()
                    B[word: ]

                TESTS:

                Try to check that :trac:`5843` Heisenbug is fixed::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: B = AlgebrasWithBasis(QQ).example(('a', 'c'))
                    sage: A == B
                    False
                    sage: Aone = A.one_from_one_basis
                    sage: Bone = B.one_from_one_basis
                    sage: Aone is Bone
                    False

               Even if called in the wrong order, they should returns their
               respective one::

                    sage: Bone().parent() is B
                    True
                    sage: Aone().parent() is A
                    True
                """
                return self.monomial(self.one_basis()) #.

            @lazy_attribute
            def one(self):
                r"""
                Return the multiplicative unit element.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word: 
                    sage: A.one()
                    B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.one_from_one_basis

            @lazy_attribute
            def from_base_ring(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.from_base_ring_from_one_basis

            def from_base_ring_from_one_basis(self, r):
                """
                Implement the canonical embedding from the ground ring.

                INPUT:

                - ``r`` -- an element of the coefficient ring

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring_from_one_basis(3)
                    3*B[word: ]
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                    sage: A(3)
                    3*B[word: ]
                """
                return self.term(self.one_basis(), r)
