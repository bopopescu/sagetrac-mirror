r"""
Complete Discrete Valuation Rings (CDVR) and Fields (CDVF)
"""
from __future__ import absolute_import
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method

from sage.categories.category_singleton import Category_singleton
from .discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
#from sage.misc.cachefunc import cached_method


class CompleteDiscreteValuationRings(Category_singleton):
    """
    The category of complete discrete valuation rings

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationRings()
        True
        sage: QQ in CompleteDiscreteValuationRings()
        False
        sage: QQ[['u']] in CompleteDiscreteValuationRings()
        True
        sage: Qp(7) in CompleteDiscreteValuationRings()
        False
        sage: TestSuite(CompleteDiscreteValuationRings()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            An error is raised when the input is indistinguishable from 0::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                Traceback (most recent call last):
                ...
                ValueError: Cannot determine the denominator of an element indistinguishable from 0
            """
            return self.parent()(1)

        def _matrix_echelonize(self, M, transformation=True, secure=False):
            """
            Row-echelonize this matrix

            INPUT:

            - ``transformation`` -- a boolean (default: True)
              Indicates whether the transformation matrix is returned

            OUTPUT:

            The position of the pivots and the transformation matrix
            if asked for.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: M.echelon_form()  # indirect doctest
                [ ...1  ...1]
                [    0 ...10]

                sage: H,L = M.echelon_form(transformation=True)  # indirect doctest
                sage: H
                [ ...1  ...1]
                [    0 ...10]
                sage: L
                [        ...1 ...444444444]
                [...444444444         ...2]
                sage: L*M == H
                True

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: M.echelon_form()  # indirect doctest
                [ ...1  ...1]
                [    0 ...10]
                [    0     0]

            An error is raised if the precision on the entries is
            not enough to determine the echelon form::

                sage: M = matrix(A, 2, 2, [A(0,5), 1, 5^8, 1])
                sage: M.echelon_form()  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to echelonize

            TESTS::

            We check that it works over various rings::

                sage: from sage.rings.padics.precision_error import PrecisionError
                sage: ring1 = ZpCA(5,15)
                sage: ring2 = Zq(5^3,names='a')
                sage: ring3 = Zp(5).extension(x^2-5, names='pi')
                sage: ring4 = PowerSeriesRing(GF(5), name='t')
                sage: for A in [ ring1, ring2, ring3, ring4 ]:
                ....:     for _ in range(10):
                ....:         M = random_matrix(A,4)
                ....:         try:
                ....:             H, L = M.echelon_form(transformation=True)
                ....:         except PrecisionError:
                ....:             continue
                ....:         if L*M != H: raise RuntimeError
            """
            from sage.matrix.matrix_cdv_dense import echelonize
            return echelonize(M, transformation, secure=secure)


        @abstract_method
        def lift_to_precision(self, absprec=None):
            """
            Return another element of the same parent with absolute precision
            at least ``absprec``, congruent to this element modulo the
            precision of this element.

            INPUT:

            - ``absprec`` -- an integer or ``None`` (default: ``None``), the
              absolute precision of the result. If ``None``, lifts to the maximum
              precision allowed.

            .. NOTE::

                If setting ``absprec`` that high would violate the precision cap,
                raises a precision error.  Note that the new digits will not
                necessarily be zero.

            EXAMPLES::

                sage: R = ZpCA(17)
                sage: R(-1,2).lift_to_precision(10)
                16 + 16*17 + O(17^10)
                sage: R(1,15).lift_to_precision(10)
                1 + O(17^15)
                sage: R(1,15).lift_to_precision(30)
                Traceback (most recent call last):
                ...
                PrecisionError: Precision higher than allowed by the precision cap.
                sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
                True

                sage: R = Zp(5); c = R(17,3); c.lift_to_precision(8)
                2 + 3*5 + O(5^8)
                sage: c.lift_to_precision().precision_relative() == R.precision_cap()
                True

            """

class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationFields()
        False
        sage: QQ in CompleteDiscreteValuationFields()
        False
        sage: LaurentSeriesRing(QQ,'u') in CompleteDiscreteValuationFields()
        True
        sage: Qp(7) in CompleteDiscreteValuationFields()
        True
        sage: TestSuite(CompleteDiscreteValuationFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]

    class ParentMethods:
        @abstract_method
        def integer_ring(self):
            """
            Return the integer ring of this CDVF

            EXAMPLES::

                sage: K = Qp(5)
                sage: K.integer_ring()
                5-adic Ring with capped relative precision 20

                sage: L.<t> = LaurentSeriesRing(GF(5))
                sage: L.integer_ring()
                Power Series Ring in t over Finite Field of size 5
            """

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            An error is raised when the input is indistinguishable from 0::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                Traceback (most recent call last):
                ...
                ValueError: Cannot determine the denominator of an element indistinguishable from 0
            """
            if self == 0:
                raise ValueError("Cannot determine the denominator of an element indistinguishable from 0")
            val = self.valuation()
            R = self.parent().integer_ring()
            if val >= 0:
                return R(1)
            else:
                return R(1) << (-val)
