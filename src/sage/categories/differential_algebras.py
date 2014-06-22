"""
Differential algebras

AUTHORS:

- Miguel Marco, John Palmieri, Travis Scrimshaw (2014-06-21): Initial version
"""
#*****************************************************************************
#  Copyright (C) 2014 Miguel Marco, John Palmieri, Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.category_types import Category_over_base_ring
from sage.categories.algebras import Algebras
from sage.categories.category_types import ChainComplexes
from sage.categories.tensor import TensorProductsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute

class DifferentialAlgebrasCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = GradedAlgebras(QQ)
            sage: C
            Category of graded algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of graded modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Graded().base_ring()
            Rational Field
            sage: GradedHopfAlgebrasWithBasis(QQ).base_ring()
            Rational Field
        """
        super(DifferentialAlgebrasCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Differential"

    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT:

        A *tuple* whose first entry is the base category class.

        .. WARNING::

            This is only used for graded categories that are not
            implemented as nested classes, and won't work otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: DifferentialAlgebras._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: DifferentialAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

                sage: F = DifferentialAlgebras
                sage: F._foo = F._base_category_class[0]
                sage: F._foo
                Traceback (most recent call last):
                ...
                AssertionError: base category class for <...Algebras'> mismatch;
                expected <...Algebras'>, got <...DifferentialAlgebras'>
        """
        module_name = cls.__module__.replace("differential_","")
        import sys
        name   = cls.__name__.replace("Differential","")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category, *args):
        """
        Magic support for putting Differential categories in their own file.

        EXAMPLES::

            sage: from sage.categoires.differential_algebras import DifferentialAlgebras
            sage: DifferentialAlgebras(ZZ)   # indirect doctest
            Category of graded modules over Integer Ring
            sage: Algebras(ZZ).Differential()
            Category of graded modules over Integer Ring
            sage: DifferentialAlgebras(ZZ) is Algebras(ZZ).Differential()
            True

        .. TODO::

            Generalize this support for all other functorial
            constructions if at some point we have a category ``Blah`` for
            which we want to implement the construction ``Blah.Foo`` in a
            separate file like we do for e.g. :class:`DifferentialAlgebras`,
            :class:`GradedAlgebras`, ...

        .. SEEALSO:: :meth:`_base_category_class`
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super(DifferentialAlgebrasCategory, cls).__classcall__(cls, category, *args)
        return base_category_class(category, *args).Differential()

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Algebra(QQ).Differnetial()  # indirect doctest
            Category of differential algebras over Rational Field
        """
        return "differential {}".format(self.base_category()._repr_object_names())

class DifferentialAlgebras(DifferentialAlgebrasCategory):
    """
    The category of differential algebras.

    EXAMPLES::

        sage: Algebras(QQ).Differential()
        Category of differential algebras with basis over Rational Field
        sage: Algebras(QQ).Differential().super_categories()
        [Category of algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(Algebras((QQ).Differential()).run()
    """
    def extra_super_categories(self):
        r"""
        Return the :class:`ChainComplexes` category.

        This method specifies that a differential algebra with
        a basis is a chain complex.

        .. SEEALSO::

            The :ref:`axioms-deduction-rules` section in the
            documentation of axioms

        EXAMPLES::

            sage: C = Algebras(QQ).Differential()
            sage: C.extra_super_categories()
            (Category of chain complexes over Rational Field,)
        """
        return [ChainComplexes(self.base_ring())]

    def example(self):
        """
        An example of differential algebra.

        EXAMPLES::

            sage: from sage.categories.differential_algebras import DifferentialAlgebras
            sage: DifferentialAlgebras(QQ).example()
            Free commutative differential algebra over Rational Field on generators in degrees 1, 1, 1, 1, 1
        """
        from sage.algebras.differential_algebras.commutative_dga import CommutativeDGA
        return CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, (((1,1,0,0,0), 1),), (((0,1,1,0,0), 1),)))

    class ParentMethods:
        def homology_algebra(self, n=None):
            """
            The homology of self as an algebra, up to but not including
            total degree ``n``

            INPUT:

            - ``n`` - integer or ``None``. If ``None``, use one more than
              the maximum of the degrees of the generators, relations,
              and images of differentials

            Return as a ``CommutativeDGA``.
            """

            # given a collection of homology generators in degrees less
            # than d, look at degree d. find a set of expected linearly
            # independent products of homology generators in that
            # degree. lift and embed into the cycles, then project to
            # homology. is the set of full rank, or are there
            # dependencies? find a basis for the subvector space of
            # homology spanned by these products, find dependencies among
            # the remaining products to find relations, and then extend to
            # a basis of the rest of homology to find new generators.

            # Generators: list of generators.

            # Relations: list of relations. Note that these are described
            # in terms of the homology generators, not in terms of the
            # generators for the original DGA.
            from sage.misc.misc_c import prod

            if n is None:
                try:
                    diffs = [self.differential(x) for x in self.gens()]
                    n = 1 + max(self._degrees + self._relations
                                + tuple([x.degree() for x in diffs if x != 0]))
                except ValueError:
                    # All lists are empty: no generators, etc.
                    n = 1

            p = self.base_ring().characteristic()
            G = self.grading_group()
            gens = []
            relns = []
            test = CommutativeDGA(self.base_ring(), degrees=())

            # in characteristic p, pth powers of algebra generators are
            # automatically homology generators.
            #
            # don't do this automatically? if an algebra generator is also
            # a homology generator, we need to then remove its pth power
            # from the list of generators, or add in a relation to
            # compensate. it may be easier just to add the pth powers as
            # they arise naturally.
            #
            # for now, leave this here to test things.
            if p > 0:
                gens = [x**p for x in self.gens() if x**p != 0]
                test = CommutativeDGA(self.base_ring(),
                                      degrees=[x.degree() for x in gens])

            update = False
            for d in group_range(1, n, group=G):
                deg = G(d)
                basis_exps = test.basis(deg)
                basis = [prod(g**e for (g,e) in zip(gens, exp))
                         for exp in basis_exps]
                H_raw = self._raw_homology(deg)
                H = self.homology(deg)

                print d, basis, H

                # under some circumstances, update either "gens" or
                # "relns", and then update "test"
                if False:
                    gens.append(blah)
                    update = True
                if False:
                    relns.append(blah)
                    update = True
                if update:
                    test = CommutativeDGA(self.base_ring(),
                                          degrees=[x.degree() for x in gens],
                                          relations=relns)
            pass

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of differential algebras constructed as Cartesian
        products of differential algebras.

        This construction gives the direct product of differential algebras.
        See discussion on:

        - http://en.wikipedia.org/wiki/Direct_product
        """
        def extra_super_categories(self):
            """
            A cartesian product of differential algebras is endowed with a
            natural algebra structure and differential.

            EXAMPLES::

                sage: C = Algebras(QQ).Differential().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of Cartesian products of commutative additive groups,
                 Category of Cartesian products of distributive magmas and additive magmas,
                 Category of Cartesian products of semigroups,
                 Category of Cartesian products of unital magmas,
                 Category of algebras over Rational Field]
            """
            return [self.base_category()]

    class TensorProducts(TensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).Differential().TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).Differential().TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of differential algebras is a
            differential algebra.
            """
            return [self.base_category()]

