r"""
Free Fermions Super Lie Conformal Algebra.

Given an `R`-module `M` with a skew-symmetric, bilinear pairing
`\langle\cdot, \cdot\rangle: M\otimes_R M \rightarrow R`. The
*Free Fermions* super Lie conformal algebra associated to this datum is
the free `R[T]`-super module generated by `\Pi M` (a purely odd copy
of `M`) plus a central vector `K` satisfying `TK=0`. The remaining
`\lambda`-brackets are given by:

.. MATH::

    [v_\lambda w] = \langle v,w \rangle K,

where `v,w \in M`.

This is an H-graded Lie conformal algebra where every generator
`v \in M` has degree `1/2`.


AUTHORS:

- Reimundo Heluani (06-03-2020): Initial implementation.
"""
#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra

class FreeFermionsLieConformalAlgebra(GradedLieConformalAlgebra):

    def __init__(self,R,gram_matrix=None,ngens=None,names=None,
                 index_set=None):
        r"""
        The Free Fermions Super Lie conformal algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``gram_matrix``: a matrix (default: ``[1]``); a symmetric
          square matrix with coefficients in ``R``.
        - ``ngens``: a positive Integer (default ``1``); the number of
          non-central generators of this Lie conformal algebra.

        OUTPUT:

        The Free Fermions Lie conformal algebra with generators
         `\psi_i`, `i=1,...,n` and `\lambda`-brackets

         .. MATH::

            [{\psi_i}_{\lambda} \psi_j] = M_{ij} K,

        where `n` is the number of generators ``ngens`` and `M` is the
        ``gram_matrix``. This super Lie conformal
        algebra is `H`-graded where every generator has degree `1/2`.

        EXAMPLES::

            sage: R = FreeFermionsLieConformalAlgebra(QQbar); R
            The free Fermions super Lie conformal algebra with generators (psi, K) over Algebraic Field.
            sage: R.inject_variables()
            Defining psi, K
            sage: psi.bracket(psi)
            {0: K}

            sage: R = FreeFermionsLieConformalAlgebra(QQbar,gram_matrix=Matrix([[0,1],[1,0]])); R
            The free Fermions super Lie conformal algebra with generators (psi_0, psi_1, K) over Algebraic Field.
            sage: R.inject_variables()
            Defining psi_0, psi_1, K
            sage: psi_0.bracket(psi_1)
            {0: K}
            sage: psi_0.degree()
            1/2
            sage: R.category()
            Category of finitely generated super H-graded Lie conformal algebras with basis over Algebraic Field
        """
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.special import identity_matrix
        if (gram_matrix is not None):
            if ngens is None:
                ngens = gram_matrix.dimensions()[0]
            try:
                assert (gram_matrix in MatrixSpace(R,ngens,ngens))
            except AssertionError:
                raise ValueError("The gram_matrix should be a symmetric " +
                    "{0} x {0} matrix, got {1}".format(ngens,gram_matrix))
            if not gram_matrix.is_symmetric():
                raise ValueError("The gram_matrix should be a symmetric " +
                    "{0} x {0} matrix, got {1}".format(ngens,gram_matrix))
        else:
            if ngens is None:
                ngens = 1;
            gram_matrix = identity_matrix(R,ngens,ngens)

        if (names is None) and (index_set is None):
            if ngens==1:
                names = 'psi'
            else:
                names = 'psi_'
            self._latex_names = tuple(r'\psi_{%d}' % i \
                                      for i in range (ngens)) + ('K',)

        from sage.structure.indexed_generators import \
                                                standardize_names_index_set
        names,index_set = standardize_names_index_set(names=names,
                                                      index_set=index_set,
                                                      ngens=ngens)
        fermiondict = { (i,j): {0: {('K',0): gram_matrix[index_set.rank(i),
                    index_set.rank(j)]}} for i in index_set for j in index_set}

        from sage.rings.rational_field import QQ
        weights = (QQ(1/2),)*ngens
        parity = (1,)*ngens
        GradedLieConformalAlgebra.__init__(self,R,fermiondict,names=names,
                                           index_set=index_set,weights=weights,
                                           parity=parity,
                                           central_elements=('K',))

        self._gram_matrix = gram_matrix

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: lie_conformal_algebras.FreeFermions(QQ)
            The free Fermions super Lie conformal algebra with generators (psi, K) over Rational Field
        """
        return "The free Fermions super Lie conformal algebra "\
                    "with generators {} over {}".format(self.gens(),
                                                         self.base_ring())

    def gram_matrix(self):
        r"""
        The Gram matrix that specifies the `\lambda`-brackets of the
        generators.

        EXAMPLES::

            sage: R = FreeFermionsLieConformalAlgebra(QQ,ngens=2);
            sage: R.gram_matrix()
            [1 0]
            [0 1]
        """
        return self._gram_matrix


