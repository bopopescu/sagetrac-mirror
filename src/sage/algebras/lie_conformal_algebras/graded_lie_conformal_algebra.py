"""
Graded Lie Conformal Algebras.

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.infinity import Infinity
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from .lie_conformal_algebra_element import GradedLCAElement
from .lie_conformal_algebra_with_structure_coefs import \
                                LieConformalAlgebraWithStructureCoefficients

class GradedLieConformalAlgebra(LieConformalAlgebraWithStructureCoefficients):
    def __init__(self, R, s_coeff, index_set=None, central_elements=None,
                 category=None, prefix=None, names=None, latex_names=None,
                 parity=None, weights=None, **kwds):
        r"""
        An H-Graded Lie conformal algebra.

        INPUT:

        - ``R`` -- a commutative ring (default: ``None``); the base
          ring of this Lie conformal algebra. Behaviour is undefined if
          it is not a field of characteristic zero

        - ``s_coeff`` -- a dictionary (default: ``None``); as in the
           input of :class:`LieConformalAlgebraStructureCoefficients`

        - ``names`` -- tuple of ``str`` (default: ``None``); as in the
          input of
          :class:`LieConformalAlgebraStructureCoefficients`

        - ``central_elements`` -- tuple of ``str`` (default: ``None``);
          as in the input of
          :class:`LieConformalAlgebraStructureCoefficients`

        - ``index_set`` -- enumerated set (default: ``None``); as in the
          input of
          :class:`LieConformalAlgebraStructureCoefficients`

        - ``weights`` -- tuple of non-negative rational numbers
          (default: tuple of ``1``); a list of degrees for this Lie
          conformal algebra.
          This tuple needs to have the same cardinality as
          ``index_set`` or ``names``. Central elements are assumed
          to have weight ``0``.

        - ``category`` The category that this Lie conformal algebra
           belongs to.

        - ``parity`` -- tuple of ``0`` or ``1`` (Default: tuple of
          ``0``); a tuple specifying the parity of each non-central
          generator.
        """
        category = LieConformalAlgebras(R).Graded().WithBasis()\
                   .FinitelyGenerated().or_subcategory(category)
        element_class = GradedLCAElement
        LieConformalAlgebraWithStructureCoefficients.__init__(self,R,
            s_coeff,index_set=index_set,central_elements=central_elements,
            category=category, element_class=element_class, prefix=prefix,
            names=names, latex_names=latex_names, parity=parity, **kwds)

        if weights is None:
            weights = (1,)* (len(self._generators) -
                             len(self.central_elements()))
        if len (weights) != (len(self._generators) -
                                len(self.central_elements())):
            raise ValueError("weights and (non-central) generator lists "\
                             "must be of same length")
        self._weights = weights

