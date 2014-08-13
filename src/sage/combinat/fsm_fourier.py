"""
Fourier Coefficients

::

    sage: import sys
    sage: sys.path.append(".")

"""

from sage.combinat.finite_state_machine import Transducer
from sage.misc.cachefunc import cached_method
class FSMFourier(Transducer):
    """
    Fourier coefficients for the sum of output of transducers.
    """

    @cached_method
    def _fourier_coefficient_data_(self):
        r"""
        Return the common data needed for the computation of all
        Fourier coefficients of the periodic fluctuation of the sum of
        output.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`namedtuple` consisting of:

        - ``c`` -- number of final components.

        - ``periods`` -- list of periods of the final components.

        - ``period`` -- least common multiple of the periods.

        - ``T`` -- eigenvector matrix.

        - ``w`` -- list of lists of vectors `\mathbf{w}_{jk}`.

        - ``coefficient_lambda`` -- list of coefficients `\lambda_j`.

        EXAMPLES:

        -   Binary sum of digits::

            sage: function('f')
            f
            sage: var('n')
            n
            sage: from fsm_fourier import FSMFourier
            sage: T = transducers.Recursion([
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            FourierCoefficientData(c=1, periods=[1], period=1, T=[1],
            w=[[(1)]], coefficient_lambda=[1])

        -   NAF::

            sage: T = transducers.Recursion([
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            (
                       [1/3   1   0]
                       [1/3   0   1]
            1, [1], 1, [1/3  -1  -1], [[(1/3, 1/3, 1/3)]], [1]
            )

        -   Abelian complexity of the paperfolding sequence::

            sage: T = transducers.Recursion([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n+2) == f(2*n+1)+1,
            ....:     f(16*n+1) == f(8*n+1),
            ....:     f(16*n+5) == f(4*n+1)+2,
            ....:     f(16*n+11) == f(4*n+3)+2,
            ....:     f(16*n+15) == f(2*n+2)+1,
            ....:     f(1) == 2, f(0) == 0]
            ....:     + [f(16*n+jj) == f(2*n+1)+2 for jj in [3,7,9,13]],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            (
                       [1/10    1    0    0    0    0    0    0    0    0]
                       [1/10    0    1    0    0    0    0    0    0    0]
                       [1/10    0    0    1    0    0    0    0    0    0]
                       [1/10    0    0    0    1    0    0    0    0    0]
                       [1/10    0    0    0    0    1    0    0    0    0]
                       [1/10    0    0    0    0    0    1    0    0    0]
                       [1/10    0    0    0    0    0    0    1    0    0]
                       [1/10    0    0    0    0    0    0    0    1    0]
                       [1/10    0    0    0    0    0    0    0    0    1]
            1, [1], 1, [1/10    0    0   -3   -2   -2   -2   -1   -1   -1],
            [[(0, 0, 3/13, 2/13, 2/13, 2/13, 1/13, 1/13, 1/13, 1/13)]],
            [1]
            )

        -   Artificial example, one-periodic, 2 states::

            sage: T = transducers.Recursion([
            ....:     f(4*n) == f(2*n)+0,
            ....:     f(4*n+2) == f(n)+1,
            ....:     f(2*n+1) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            (
                       [1/2   1]
            1, [1], 1, [1/2  -1],
            [[(1/2, 1/2)]],
            [1]
            )

        -   Artificial example, period 3::

            sage: T = transducers.Recursion([
            ....:     f(8*n) == f(4*n+3)+3,
            ....:     f(8*n+4) == f(4*n+3)+1,
            ....:     f(8*n+2) == f(4*n+3)+2,
            ....:     f(8*n+6) == f(4*n+3)-1,
            ....:     f(8*n+1) == f(4*n)+5,
            ....:     f(8*n+5) == f(4*n+2)+1,
            ....:     f(8*n+3) == f(4*n+1)+2,
            ....:     f(8*n+7) == f(4*n+1),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            (
                          [         1/7            1            1            1            0            0            0]
                          [         1/7      4*zeta3 -4*zeta3 - 4            0            1            0            0]
                          [         1/7     -2*zeta3  2*zeta3 + 2            0            0            1            0]
                          [         1/7 -4*zeta3 - 4      4*zeta3            0            0            0            1]
                          [         1/7 -4*zeta3 - 4      4*zeta3            0            0            0           -1]
                          [         1/7      4*zeta3 -4*zeta3 - 4            0            0            0            0]
               1, [3], 3, [         1/7            4            4            0            0            0            0],
               [[(0, 0, 0, 1/6, 1/6, 1/3, 1/3),
                 (0, 0, 0, 1/24*zeta3, 1/24*zeta3, -1/12*zeta3 - 1/12, 1/12),
                 (0, 0, 0, -1/24*zeta3 - 1/24, -1/24*zeta3 - 1/24, 1/12*zeta3, 1/12)]],
               [1]
            )

        -   Artificial example, period 2, vanishing w-vector::

            sage: T = transducers.Recursion([
            ....:     f(4*n) == f(2*n+1)+1,
            ....:     f(4*n+1) == f(2*n)+2,
            ....:     f(4*n+2) == f(2*n+1)+3,
            ....:     f(4*n+3) == f(2*n)-1,
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: FSMFourier(T)._fourier_coefficient_data_()
            (
                       [1/3   0   1]
                       [1/3   1   0]
            1, [2], 2, [1/3  -1   0],
            [[(0, 1/2, 1/2), (0, 0, 0)]],
            [1]
            )

        -   Artificial example with two final components of periods `2`
            and `3`, respectively::

            sage: T = FSMFourier([(0, 1, 0, 1), (1, 2, 0, 1),
            ....:     (2, 1, 0, 2), (1, 2, 1, 0), (2, 1, 1, 2),
            ....:     (0, -1, 1, 1), (-1, -2, 1, 1), (-2, -3, 1, 1),
            ....:     (-3, -1, 1, 1), (-1, -2, 0, 2), (-2, -3, 0, 1),
            ....:     (-3, -1, 0, 3)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2, -3, -2, -1])
            sage: T._fourier_coefficient_data_()
            (
                          [        1/7           1           1         1/5           1           1]
                          [          0           0           0         2/5          -2           0]
                          [          0           0           0         2/5           2           0]
                          [        2/7           2           2           0           0           0]
                          [        2/7    -2*zeta6 2*zeta6 - 2           0           0           0]
            2, [3, 2], 6, [        2/7 2*zeta6 - 2    -2*zeta6           0           0           0],
            [[(0, 0, 0, 1/6, 1/6, 1/6),
              (0, 0, 0, 1/6, 1/6*zeta6 - 1/6, -1/6*zeta6),
              (0, 0, 0, 1/6, -1/6*zeta6, 1/6*zeta6 - 1/6)],
             [(0, 1/4, 1/4, 0, 0, 0), (0, -1/4, 1/4, 0, 0, 0)]],
            [1/2, 1/2]
            )
        """
        import collections
        import itertools
        import operator

        from sage.matrix.constructor import matrix
        from sage.modules.free_module import VectorSpace
        from sage.modules.free_module_element import vector
        from sage.rings.arith import lcm
        from sage.rings.integer_ring import ZZ
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.structure.sage_object import SageObject

        FourierCoefficientData = collections.namedtuple(
            "FourierCoefficientData",
            ["c", "periods", "period", "T", "w", "coefficient_lambda"])

        positions = dict((state.label(), j)
                         for j, state in enumerate(self.iter_states()))
        q = len(self.input_alphabet)

        class FCComponent(SageObject):
            """Hold a final component and associated data."""
            def __init__(self, fsm):
                self.fsm = fsm
                self.period = fsm.graph().period()
                self.n_states = len(self.fsm.states())

            def eigenvectors(self, M, components, common_period):
                nrows = sum(c.n_states for c in components if c != self)
                mask = matrix(
                    nrows, M.ncols(),
                    [standard_basis[positions[state.label()]]
                     for other in components
                     if other != self
                     for state in other.fsm.iter_states()])

                def eigenvector(j):
                    eigenvalue = q * alpha**(j * common_period / self.period)
                    S = matrix.block(
                        [[M - eigenvalue*matrix.identity(M.nrows())],
                         [mask]],
                        subdivide=False)
                    kernel = S.right_kernel_matrix()
                    assert kernel.nrows() == 1
                    if j == 0:
                        #normalize for positive eigenvector
                        return kernel.row(0) / sum(kernel.row(0))
                    else:
                        return kernel.row(0)

                return [eigenvector(j) for j in range(self.period)]

            @cached_method()
            def right_eigenvectors(self):
                return self.eigenvectors(M, components, period)

            @cached_method()
            def left_eigenvectors(self):
                left_eigenvectors = self.eigenvectors(M.transpose(), components, period)
                return [w/(v*w) for v, w
                        in itertools.izip(self.right_eigenvectors(),
                                          left_eigenvectors)]

            @cached_method()
            def vectors_w(self):
                return [(initial_vector*v)*w for v, w
                        in itertools.izip(self.right_eigenvectors(),
                                          self.left_eigenvectors())]

            @cached_method()
            def coefficient_lambda(self):
                ones = vector(1 for _ in range(M.nrows()))
                products = [w*ones for w in self.vectors_w()]
                assert all(e.is_zero() for e in products[1:])
                return products[0]


        components = [FCComponent(c) for c in self.final_components()]
        period = lcm([c.period for c in components])
        field = CyclotomicField(period)
        alpha = field.gen()
        M = self.adjacency_matrix(entry=lambda t: 1)
        standard_basis = VectorSpace(field, M.nrows()).basis()

        if len(self.initial_states()) != 1:
            raise NotImplementedError(
                "Transducer does not have a unique initial state.")
        initial_vector = standard_basis[positions[
                self.initial_states()[0].label()]]

        right_eigenvectors = list(itertools.chain(
                *(c.right_eigenvectors()
                  for c in components)))

        left_eigenvectors = list(itertools.chain(
                *(c.left_eigenvectors()
                  for c in components)))

        annihilated_by_left = matrix(left_eigenvectors).\
            right_kernel_matrix().transpose()

        T = matrix.block([[matrix.column(right_eigenvectors),
                           annihilated_by_left]],
                         subdivide=False)

        assert T.is_square()
        assert T.nrows() == M.nrows()
        assert T.is_invertible()

        check = T.inverse() * M * T
        eigenvalues = [q * alpha**(j * period/c.period)
                       for c in components
                       for j in range(c.period)]
        check_dont_care = check.submatrix(len(eigenvalues),
                                          len(eigenvalues))
        assert (matrix.block(
                [[matrix.diagonal(eigenvalues), ZZ(0)],
                 [ZZ(0), check_dont_care]],
                     subdivide=False) - check).is_zero()

        assert (T.inverse().submatrix(nrows=len(left_eigenvectors))
                - matrix(left_eigenvectors)).is_zero()

        return FourierCoefficientData(
            c=len(components),
            periods=[c.period for c in components],
            period=period,
            T=T,
            w=[c.vectors_w() for c in components],
            coefficient_lambda=[c.coefficient_lambda() for c in components])
