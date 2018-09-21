r"""
Finite `\ZZ`-modules with with bilinear and quadratic forms.

AUTHORS:

- Simon Brandhorst (2017-09): First created
"""

#*****************************************************************************
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.fg_pid.fgp_element import DEBUG, FGP_Element
from sage.modules.free_quadratic_module import FreeQuadraticModule
from sage.arith.misc import gcd
from sage.rings.all import ZZ, Zp, QQ, IntegerModRing
from sage.groups.additive_abelian.qmodnz import QmodnZ
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.arith.misc import legendre_symbol
from sage.rings.finite_rings.integer_mod import mod

def TorsionQuadraticForm(q):
    r"""
    Create a torsion quadratic form module from a rational matrix.

    The resulting quadratic form takes values in `\QQ / \ZZ`
    or `\QQ / 2 \ZZ` (depending on ``q``).
    If it takes values modulo `2`, then it is non-degenerate.
    In any case the bilinear form is non-degenerate.

    INPUT:

    - ``q`` -- a symmetric rational matrix

    EXAMPLES::

        sage: q1 = Matrix(QQ,2,[1,1/2,1/2,1])
        sage: TorsionQuadraticForm(q1)
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [  1 1/2]
        [1/2   1]

    In the following example the quadratic form is degenerate.
    But the bilinear form is still non-degenerate::

        sage: q2 = diagonal_matrix(QQ,[1/4,1/3])
        sage: TorsionQuadraticForm(q2)
        Finite quadratic module over Integer Ring with invariants (12,)
        Gram matrix of the quadratic form with values in Q/Z:
        [7/12]
    """
    q = matrix(QQ, q)
    if q.nrows() != q.ncols():
        raise ValueError("the input must be a square matrix")
    if q != q.transpose():
        raise ValueError("the input must be a symmetric matrix")

    Q, d = q._clear_denom()
    S, U, V = Q.smith_form()
    D = U * q * V
    Q = FreeQuadraticModule(ZZ, q.ncols(), inner_product_matrix=d**2 * q)
    denoms = [D[i,i].denominator() for i in range(D.ncols())]
    rels = Q.span(diagonal_matrix(ZZ, denoms) * U)
    return TorsionQuadraticModule((1/d)*Q, (1/d)*rels, modulus=1)

class TorsionQuadraticModuleElement(FGP_Element):
    r"""
    An element of a torsion quadratic module.

    INPUT:

    - ``parent`` -- parent

    - ``x`` -- element of ``parent.V()``

    - ``check`` -- bool (default: ``True``)

    TESTS::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
        sage: loads(dumps(T)) == T
        True
        sage: t = T.gen(0)
        sage: loads(dumps(t)) == t
        True

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
        sage: b = V.basis()
        sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
        sage: Q = TorsionQuadraticModule(V, W)
        sage: x = Q(b[0] - b[1])
        sage: TestSuite(x).run()
        """

    def _mul_(self, other):
        r"""
        Compute the inner product of two elements.

        OUTPUT:

        - an element of `\QQ / m\ZZ` with `m\ZZ = (V, W)`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = (1/2)*ZZ^2; W = ZZ^2
            sage: T = TorsionQuadraticModule(V, W)
            sage: g = T.gens()
            sage: x = g[0]
            sage: y = g[0] + g[1]
            sage: x
            (1, 0)
            sage: x*y
            1/4

        The inner product has further aliases::

            sage: x.inner_product(y)
            1/4
            sage: x.b(y)
            1/4
        """
        value_module = self.parent().value_module()
        return value_module( self.lift().inner_product(other.lift()) )

    inner_product = _mul_
    b = _mul_

    def quadratic_product(self):
        r"""
        Compute the quadratic_product of ``self``.

        OUTPUT:

        - an element of `\QQ / n\ZZ` where `n\ZZ = 2(V,W) +
          \ZZ \{ (w,w) | w \in W \}`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: W = FreeQuadraticModule(ZZ, 2, 2*matrix.identity(2))
            sage: V = (1/2) * W
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.gen(0)
            sage: x
            (1, 0)
            sage: x.quadratic_product()
            1/2
            sage: x.quadratic_product().parent()
            Q/2Z
            sage: x*x
            1/2
            sage: (x*x).parent()
            Q/Z
        """
        value_module_qf = self.parent().value_module_qf()
        lift = self.lift()
        return value_module_qf(lift.inner_product(lift))

    q = quadratic_product

class TorsionQuadraticModule(FGP_Module_class):
    r"""
    Finite quotients with a bilinear and a quadratic form.

    Let `V` be a symmetric FreeQuadraticModule and `W \subseteq V` a
    submodule of the same rank as `V`. The quotient `V / W` is a torsion
    quadratic module. It inherits a bilinear form `b` and a quadratic
    form `q`.

    `b: V \times V \to \QQ / m\ZZ`, where  `m\ZZ = (V,W)`
    and `b(x,y) = (x,y) + m\ZZ`

    `q: V \to \QQ / n\ZZ`, where `n\ZZ = 2(V,W) + \ZZ \{ (w,w) | w \in W \}`

    INPUT:

    - ``V`` -- a :class:`FreeModule` with a symmetric inner product matrix

    - ``W`` -- a submodule of ``V`` of the same rank as ``V``

    - ``check`` -- bool (default: ``True``)

    - ``modulus`` -- a rational number dividing `m` (default: `m`);
      the inner product `b` is defined in `\QQ /` ``modulus`` `\ZZ`

    - ``modulus_qf`` -- a rational number dividing `n` (default: `n`);
      the quadratic form `q` is defined in `\QQ /` ``modulus_qf`` `\ZZ`

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = FreeModule(ZZ, 3)
        sage: T = TorsionQuadraticModule(V, 5*V)
        sage: T
        Finite quadratic module over Integer Ring with invariants (5, 5, 5)
        Gram matrix of the quadratic form with values in Q/5Z:
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    Element = TorsionQuadraticModuleElement

    def __init__(self, V, W, gens=None, modulus=None, modulus_qf=None, check=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
            sage: TestSuite(T).run()
        """
        if check:
            if V.rank() != W.rank():
                raise ValueError("modules must be of the same rank")
            if V.base_ring() is not ZZ:
                raise NotImplementedError("only currently implemented over ZZ")
            if V.inner_product_matrix() != V.inner_product_matrix().transpose():
                raise ValueError("the cover must have a symmetric inner product")

            if gens is not None and V.span(gens) + W != V:
                raise ValueError("provided gens do not generate the quotient")

        FGP_Module_class.__init__(self, V, W, check=check)
        if gens is not None:
            self._gens_user = tuple(self(v) for v in gens)
        else:
            # this is taken care of in the .gens method
            # we do not want this at initialization
            self._gens_user = None

        # compute the modulus - this may be expensive
        if modulus is None or check:
            # The inner product of two elements `b(v1+W,v2+W)`
            # is defined `mod (V,W)`
            num = V.basis_matrix() * V.inner_product_matrix() * W.basis_matrix().T
            self._modulus = gcd(num.list())

        if modulus is not None:
            if check and self._modulus / modulus not in self.base_ring():
                raise ValueError("the modulus must divide (V, W)")
            self._modulus = modulus

        if modulus_qf is None or check:
            # The quadratic_product of an element `q(v+W)` is defined
            # `\mod 2(V,W) + ZZ\{ (w,w) | w in w\}`
            norm = gcd(self.W().gram_matrix().diagonal())
            self._modulus_qf = gcd(norm, 2 * self._modulus)

        if modulus_qf is not None:
            if check and self._modulus_qf / modulus_qf not in self.base_ring():
                raise ValueError("the modulus_qf must divide (V, W)")
            self._modulus_qf = modulus_qf

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ,3)
            sage: T = TorsionQuadraticModule(V, 5*V,modulus=1)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        return ( "Finite quadratic module over %s with invariants %s\n"
                 % (self.base_ring(),self.invariants()) +
                 "Gram matrix of the quadratic form with values in %r:\n%r"
                 % (self.value_module_qf(),self.gram_matrix_quadratic()) )

    def _module_constructor(self, V, W, check=False):
        r"""
        Construct a torsion quadratic module ``V / W``.

        INPUT:

        - ``V`` -- an module

        - ``W`` -- an submodule of ``V`` over the same base ring

        - ``check`` -- bool (default: ``False``);

          * if ``False``, then the value modulus is inherited from ``self``
          * if ``True``, it figures it out on its own. But that is expensive

        OUTPUT:

        The quotient ``V / W`` as a :class:`TorsionQuadraticModule`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
            sage: b = V.basis()
            sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
            sage: Q = TorsionQuadraticModule(V, W); Q
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
            sage: Q._module_constructor(V,W)
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
        """
        if check:
            # figuring out the modulus can be expensive
            return TorsionQuadraticModule(V, W, check=check)
        else:
            return TorsionQuadraticModule(V, W, check=check,
                                          modulus=self._modulus,
                                          modulus_qf=self._modulus_qf)

    @cached_method
    def _to_smith(self):
        r"""
        Return the tranformation matrix from generators to smith generators.

        EXAMPLES::
        """
        # transformation matrices between smith and normal_form generators
        to_smith = matrix(ZZ, [t.vector() for t in self.gens()])
        return to_smith

    @cached_method
    def _to_gens(self):
        r"""
        Return the transformation matrix from smith generators to generators.

        EXAMPLES::

            sage:
        """
        invs = self.invariants()
        if len(invs) == 0:
            return matrix([])
        R = IntegerModRing(invs[-1])
        E = matrix.identity(R, len(invs))
        # view self as a submodule of (ZZ/nZZ)^(len(invs))
        B = self._to_smith().change_ring(R)
        for k in range(B.ncols()):
            B[:, k] *= invs[-1] // invs[k]
            E[:, k] *= invs[-1] // invs[k]
        # unfortunatly solve_left does not work with matrices
        # to_normal = matrix([B.solve_left(e) for e in E.rows()])
        to_normal = B.solve_left(E).change_ring(ZZ)
        return to_normal

    def all_submodules(self):
        r"""
        Return a list of all submodules of ``self``.

        WARNING:

        This method creates all submodules in memory. The number of submodules
        grows rapidly with the number of generators. For example consider a
        vector space of dimension `n` over a finite field of prime order `p`.
        The number of subspaces is (very) roughly `p^{(n^2-n)/2}`.

        EXAMPLES::

            sage: D = IntegralLattice("D4").discriminant_group()
            sage: D.all_submodules()
            [Finite quadratic module over Integer Ring with invariants ()
              Gram matrix of the quadratic form with values in Q/2Z:
              [],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2, 2)
              Gram matrix of the quadratic form with values in Q/2Z:
              [  1 1/2]
              [1/2   1]]
        """
        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        invs = self.invariants()
        # knows how to compute all subgroups
        A = AbelianGroupGap(invs)
        S = A.all_subgroups()
        # over ZZ submodules and subgroups are the same thing.
        submodules = []
        for sub in S:
            gen = [A(g).exponents() for g in sub.gens()]
            gen = [self.linear_combination_of_smith_form_gens(g) for g in gen]
            submodules.append(self.submodule(gen))
        return submodules

    def all_totally_isotropic_submodules(self, bilinear=False):
        r"""
        Return a list of all totally isotropic submodules.

        Be aware that the number of such submodules grows quickly.
        The method works by first creating all subgroups and then discarding the
        non - isotropic ones. This can be time and memory consuming.

        EXAMPLES::

            sage: q = Matrix(QQ,2,[2,1,1,2])/2
            sage: T = TorsionQuadraticForm(q)
            sage: T.all_totally_isotropic_submodules()
            [Finite quadratic module over Integer Ring with invariants ()
            Gram matrix of the quadratic form with values in Q/2Z:
            []]
            sage: T.all_totally_isotropic_submodules(bilinear=True)
            [Finite quadratic module over Integer Ring with invariants ()
            Gram matrix of the quadratic form with values in Q/2Z:
            [], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1], Finite quadratic module over Integer Ring with invariants (2,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1]]
        """
        isotropic = []
        n = self.value_module_qf().n
        for sub in self.all_submodules():
            if bilinear:
                q = sub.gram_matrix_bilinear()
            else:
                q = sub.gram_matrix_quadratic()
            if q == 0:
                isotropic.append(sub)
        return isotropic

    @cached_method
    def brown_invariant(self):
        r"""
        Return the Brown invariant of this torsion quadratic form.

        Let `(D,q)` be a torsion quadratic module with values in `\QQ / 2 \ZZ`.
        The Brown invariant `Br(D,q) \in \Zmod{8}` is defined by the equation

        .. MATH::

            \exp \left( \frac{2 \pi i }{8} Br(q)\right) =
            \frac{1}{\sqrt{D}} \sum_{x \in D} \exp(i \pi q(x)).

        The Brown invariant is additive with respect to direct sums of
        torsion quadratic modules.

        OUTPUT:

        - an element of `\Zmod{8}`

        EXAMPLES::

            sage: L = IntegralLattice("D4")
            sage: D = L.discriminant_group()
            sage: D.brown_invariant()
            4

        We require the quadratic form to be defined modulo `2 \ZZ`::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3))
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: T.brown_invariant()
            Traceback (most recent call last):
            ...
            ValueError: the torsion quadratic form must have values in QQ / 2 ZZ
        """
        if self._modulus_qf != 2:
            raise ValueError("the torsion quadratic form must have values in "
                             "QQ / 2 ZZ")
        from sage.quadratic_forms.genera.normal_form import collect_small_blocks
        brown = IntegerModRing(8).zero()
        for p in self.annihilator().gen().prime_divisors():
            q = self.primary_part(p).normal_form()
            q = q.gram_matrix_quadratic()
            L = collect_small_blocks(q)
            for qi in L:
                brown += _brown_indecomposable(qi,p)
        return brown

    @cached_method
    def gram_matrix_bilinear(self):
        r"""
        Return the gram matrix with respect to the generators.

        OUTPUT:

        A rational matrix ``G`` with ``G[i,j]`` given by the inner product
        of the `i`-th and `j`-th generator. Its entries are only well
        defined `\mod (V, W)`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ, 3, matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.gram_matrix_bilinear()
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i+1):
                G[i,j] = G[j,i] = (gens[i] * gens[j]).lift()
        return G

    @cached_method
    def gram_matrix_quadratic(self):
        r"""
        The gram matrix of the quadratic form with respect to the generators.

        OUTPUT:

        - a rational matrix ``Gq`` with ``Gq[i,j] = gens[i]*gens[j]``
          and ``G[i,i] = gens[i].q()``

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ, [[2,0,0,-1],[0,2,0,-1],[0,0,2,-1],[-1,-1,-1,2]])
            sage: D4 = FreeQuadraticModule(ZZ, 4, D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: discrForm = TorsionQuadraticModule(D4dual, D4)
            sage: discrForm.gram_matrix_quadratic()
            [  1 1/2]
            [1/2   1]
            sage: discrForm.gram_matrix_bilinear()
            [  0 1/2]
            [1/2   0]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i):
                G[i,j] = G[j,i] = (gens[i] * gens[j]).lift()
            G[i,i] = gens[i].q().lift()
        return G

    def gens(self):
        r"""
        Return generators of ``self``.

        There is no assumption on the generators except that they
        generate the module.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 3)
            sage: T = TorsionQuadraticModule(V, 5*V)
            sage: T.gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        if self._gens_user is None:
            return self.smith_form_gens()
        return self._gens_user

    def genus(self, signature_pair):
        r"""
        Return the genus defined by ``self`` and the ``signature_pair``.

        EXAMPLES::

            sage: L = IntegralLattice("D4")
            sage: D = L.discriminant_group()
            sage: Gen = D.genus(L.signature_pair())
            sage: Gen == L.genus()
            True
        """
        from sage.quadratic_forms.genera.genus import (Genus_Symbol_p_adic_ring,
                                                    GenusSymbol_global_ring,
                                                    p_adic_symbol,
                                                    is_GlobalGenus)
        from sage.misc.misc_c import prod
        s_plus = signature_pair[0]
        s_minus = signature_pair[1]
        rank = s_plus + s_minus
        D = self.cardinality()
        determinant = (-1)**s_minus * D
        local_symbols = []
        for p in (2*D).prime_divisors():
            D = self.primary_part(p)
            if len(D.invariants()) != 0:
                G_p = D.gram_matrix_quadratic().inverse()
                # get rid of denominators without changeing the local equivalence class
                G_p *= G_p.denominator()**2
                G_p = G_p.change_ring(ZZ)
                local_symbol = p_adic_symbol(G_p, p, D.invariants()[-1].valuation(p))
            else:
                local_symbol = []

            rk = rank - len(D.invariants())
            if rk > 0:
                if p == 2:
                    det = determinant.prime_to_m_part(2) % 8
                    det *= prod([di[2] for di in local_symbol])
                    det = det % 8
                    local_symbol.append([0, rk, det, 0, 0])
                else:
                    det = legendre_symbol(determinant.prime_to_m_part(p),p)
                    det = (det * prod([di[2] for di in local_symbol]))
                    local_symbol.append([0, rk, det])
            local_symbol.sort()
            local_symbol = Genus_Symbol_p_adic_ring(p, local_symbol)
            local_symbols.append(local_symbol)

        # This genus has the right discriminant group
        # but it may be empty
        genus = GenusSymbol_global_ring(signature_pair, local_symbols)

        if self.value_module_qf().n != 2:
            raise NotImplementedError(
                "Currently, this is only implemented for even genera. " +
                "Want to help us implement this for odd lattices?")

        # the symbol for p=2 and scale 1 is only well defined mod 4
        # when a jordan block of scale 1 is odd.
        # In this case the symbol is determined by the property
        # that it forms a valid global genus symbol.
        # We simply try out all possibilities untill we reach a valid symbol
        if is_GlobalGenus(genus):
            return genus

        sym2 = local_symbols[0].symbol_tuple_list()

        if sym2[0][0] == 1:
            i = 0
        else:
            i = 1

        # wrong oddity
        sym2[i][4] = (sym2[i][4] + 4) % 8
        local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
        if is_GlobalGenus(genus):
            return genus

        # wrong oddity and det
        if i == 1:
            # if there is a part of scale 0,
            # then that has to change det as well
            sym2[0][2] = (sym2[0][2] + 4) % 8
        sym2[i][2] = (sym2[i][2] + 4) % 8
        local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
        if is_GlobalGenus(genus):
            return genus

        # wrong det
        sym2[i][4] = (sym2[i][4] + 4) % 8
        local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
        if is_GlobalGenus(genus):
            return genus
        else:
            raise AssertionError("oops")

    @cached_method
    def is_degenerate(self):
        r"""
        Return if the underlying bilinear form is degenerate.

        EXAMLES::

            sage: T = TorsionQuadraticForm(matrix([1/27]))
            sage: D = T.submodule([T.gen(0)*3])
            sage: D.is_degenerate()
            True
        """
        return 1 != self.orthogonal_submodule_to(self.gens()).cardinality()

    def is_genus(self, signature_pair, even=True):
        r"""
        Return ``True`` if there is a lattice with this signature and discriminant form.

        TODO:

        - implement the same for odd lattices

        INPUT:

        - signature_pair -- a tuple of non negative integers ``(s_plus, s_minus)``
        - even -- bool (default: ``True``)

        EXAMPLES::

            sage: L = IntegralLattice("D4").direct_sum(IntegralLattice(3 * Matrix(ZZ,2,[2,1,1,2])))
            sage: D = L.discriminant_group()
            sage: D.is_genus((6,0))
            True

        Let us see if there is a lattice in the genus defined by the same discriminant form
        but with a different signature::

            sage: D.is_genus((4,2))
            False
            sage: D.is_genus((16,2))
            True
        """
        s_plus = ZZ(signature_pair[0])
        s_minus = ZZ(signature_pair[1])
        if s_plus<0 or s_minus<0:
            raise ValueError("signature invariants must be non negative")
        rank = s_plus + s_minus
        signature = s_plus - s_minus
        D = self.cardinality()
        det = (-1)**s_minus * D
        if rank < len(self.invariants()):
            return False
        if even and self._modulus_qf != 2:
            raise ValueError("the discriminant form of an even lattice has"
                                 "values modulo 2.")
        if (not even) and not (self._modulus == self._modulus_qf == 1):
            raise ValueError("the discriminant form of an odd lattice has"
                             "values modulo 1.")
        if not even:
            raise NotImplementedError("at the moment sage knows how to do this only for even genera. " +
                                      " Help us to implement this for odd genera.")
        for p in D.prime_divisors():
            # check the determinat conditions
            Q_p = self.primary_part(p)
            gram_p = Q_p.gram_matrix_quadratic()
            length_p = len(Q_p.invariants())
            u = det.prime_to_m_part(p)
            up = gram_p.det().numerator().prime_to_m_part(p)
            if p!=2 and length_p==rank:
                if legendre_symbol(u, p) != legendre_symbol(up, p):
                    return False
            if p == 2:
                if rank % 2 != length_p % 2:
                    return False
                n = (rank - length_p) / 2
                if u % 4 != (-1)**(n % 2) * up % 4:
                    return False
                if rank == length_p:
                    a = QQ(1) / QQ(2)
                    b = QQ(3) / QQ(2)
                    diag = gram_p.diagonal()
                    if not (a in diag or b in diag):
                        if u % 8 !=  up % 8:
                            return False
        if self.brown_invariant() != signature:
            return False
        return True

    def is_isomorphic_to(self, other):
        r"""
        Return if the underlying quadratic forms are isomorphic.
        """
        if self._modulus_qf != 2:
            raise NotImplementedError()
        if self.value_module_qf() != other.value_module():
            return False
        if self.invariants() != other.invariants():
            return False
        if self.is_degenerate() != other.is_degenerate():
            return False
        if len(self.invariants())==0:
            return True
        qf1 = self.normal_form().gram_matrix_quadratic()
        qf2 = other.normal_form().gram_matrix_quadratic()
        if (qf1 - qf2).denominator() != 1: # the bilinear forms agree
            return False
        if not self.is_degenerate():
            if ZZ(2).divides(self.invariants()[-1]):
                return qf1 == qf2 # check the diagonal agrees too
            return True
        if len(set(self.invariants())) == 1 and not ZZ(2).divides(self.invariants()[-1]):
            return True
        # now the normal form is not unique anymore.
        ker1 = self.orthogonal_submodule_to(self.gens())
        ker2 = other.orthogonal_submodule_to(other.gens())
        if ker1.invariants() != ker2.invariants():
            return False
        if 2 != self.invariants()[-1]:
            raise NotImplementedError(self)
        d1 = set(ker1.gram_matrix_quadratic().diagonal())
        d2 = set(ker2.gram_matrix_quadratic().diagonal())
        return d1 == d2


    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule, list, or tuple of generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 10)
            sage: T = TorsionQuadraticModule(V, 3*V)
            sage: S = T.submodule(T.gens()[:5])
            sage: O = T.orthogonal_submodule_to(S)
            sage: O
            Finite quadratic module over Integer Ring with invariants (3, 3, 3, 3, 3)
            Gram matrix of the quadratic form with values in Q/3Z:
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: O.V() + S.V() == T.V()
            True
        """
        if not isinstance(S,TorsionQuadraticModule):
            S = self.submodule(S)
        else:
            if not S.is_submodule(self):
                raise ValueError("S must be a submodule of this module")

        G = self.V().inner_product_matrix()
        T = self.V().basis_matrix()
        S = S.V().basis_matrix()
        m = self._modulus

        Y = T * G * S.transpose()
        # Elements of the ambient module which pair integrally with self.V()
        integral = Y.inverse() * T
        # Element of the ambient module which pair in mZZ with self.V()
        orthogonal = m * integral
        orthogonal = self.V().span(orthogonal.rows())
        # We have to make sure we get a submodule
        orthogonal = orthogonal.intersection(self.V())
        orthogonal = self.submodule(orthogonal.gens())
        return orthogonal

    @cached_method
    def normal_form(self, partial=False):
        r"""
        Return the normal form of this torsion quadratic module.

        Two torsion quadratic modules are isomorphic if and only if they have
        the same value modules and the same normal form.

        A torsion quadratic module `(T,q)` with values in `\QQ/n\ZZ` is
        in normal form if the rescaled quadratic module `(T, q/n)`
        with values in `\QQ/\ZZ` is in normal form.

        For the definition of normal form see [MirMor2009]_ IV Definition 4.6.
        Below are some of its properties.
        Let `p` be odd and `u` be the smallest non-square modulo `p`.
        The normal form is a diagonal matrix with diagonal entries either `p^n`
        or `u p^n`.

        If `p = 2` is even, then the normal form consists of
        1 x 1 blocks of the form

        .. MATH::

            (0), \quad 2^n(1),\quad 2^n(3),\quad 2^n(5) ,\quad 2^n(7)

        or of `2 \times 2` blocks of the form

        .. MATH::

            2^n
            \left(\begin{matrix}
                2 & 1\\
                1 & 2
            \end{matrix}\right), \quad
            2^n
            \left(\begin{matrix}
                0 & 1\\
                1 & 0
            \end{matrix}\right).

       The blocks are ordered by their valuation.

        INPUT:

        - partial - bool (default: ``False``) return only a partial normal form
          it is not unique but still useful to extract invariants

        OUTPUT:

        - a torsion quadratic module

        EXAMPLES::

            sage: L1=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,4]]))
            sage: L1.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]
            sage: L2=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: L2.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]

        We check that :trac:`24864` is fixed::

            sage: L1=IntegralLattice(matrix([[-4,0,0],[0,4,0],[0,0,-2]]))
            sage: AL1=L1.discriminant_group()
            sage: L2=IntegralLattice(matrix([[-4,0,0],[0,-4,0],[0,0,2]]))
            sage: AL2=L2.discriminant_group()
            sage: AL1.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]
            sage: AL2.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]

        Some exotic cases::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2])
            sage: D4 = FreeQuadraticModule(ZZ,4,D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: T = TorsionQuadraticModule((1/6)*D4dual,D4)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/18 5/36    0    0]
            [5/36 1/18 5/36 5/36]
            [   0 5/36 1/36 1/72]
            [   0 5/36 1/72 1/36]
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [ 1/6 1/12    0    0    0    0    0    0]
            [1/12  1/6    0    0    0    0    0    0]
            [   0    0 1/12 1/24    0    0    0    0]
            [   0    0 1/24 1/12    0    0    0    0]
            [   0    0    0    0  1/9    0    0    0]
            [   0    0    0    0    0  1/9    0    0]
            [   0    0    0    0    0    0  1/9    0]
            [   0    0    0    0    0    0    0  1/9]

        TESTS:

        A degenerate case::

            sage: T = TorsionQuadraticModule((1/6)*D4dual, D4, modulus=1/36)
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/18)Z:
            [1/36 1/72    0    0    0    0    0    0]
            [1/72 1/36    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
        """
        gens = []
        from sage.quadratic_forms.genera.normal_form import p_adic_normal_form, _normalize
        for p in self.annihilator().gen().prime_divisors():
            D_p = self.primary_part(p)
            q_p = D_p.gram_matrix_quadratic()
            q_p = q_p / D_p._modulus_qf

            # continue with the non-degenerate part
            r = q_p.rank()
            if r != q_p.ncols():
                U = q_p._clear_denom()[0].hermite_form(transformation=True)[1]
            else:
                U = q_p.parent().identity_matrix()
            kernel = U[r:,:]
            nondeg = U[:r,:]
            q_p = nondeg * q_p * nondeg.T

            # the normal form is implemented for p-adic lattices
            # so we should work with the lattice q_p --> q_p^-1
            q_p1 = q_p.inverse()
            prec = self.annihilator().gen().valuation(p) + 5
            D, U = p_adic_normal_form(q_p1, p, precision=2*prec + 5, partial=partial)
            # if we compute the inverse in the p-adics everything explodes --> go to ZZ
            U = U.change_ring(ZZ).inverse().transpose()

            # the inverse is in normal form - so to get a normal form for the original one
            # it is enough to massage each 1x1 resp. 2x2 block.
            U = U.change_ring(Zp(p, type='fixed-mod', prec=prec)).change_ring(ZZ)
            D = U * q_p * U.T * p**q_p.denominator().valuation(p)
            D = D.change_ring(Zp(p, type='fixed-mod', prec=prec))
            _, U1 = _normalize(D, normal_odd=False)
            U = U1.change_ring(ZZ) * U

            # reattach the degenerate part
            nondeg = U * nondeg
            U = nondeg.stack(kernel)

            #apply U to the generators
            n = U.ncols()
            gens_p = []
            for i in range(n):
                g = self.V().zero()
                for j in range(n):
                    g += D_p.gens()[j].lift() * U[i,j]
                gens_p.append(g)
            gens += gens_p
        return self.submodule_with_gens(gens)

    @cached_method
    def orthogonal_group(self):
        r"""
        Orthognal group of the associated torsion quadratic form.

        WARNING:

        This is usually smaller than the orthogonal group of the bilinear form.

        EXAMPLES::

            sage: L = IntegralLattice(2*matrix.identity(3))
            sage: D = L.discriminant_group()
            sage: OD = D.orthogonal_group()
            sage: OD
            Group of isometries of
            Finite quadratic module over Integer Ring with invariants (2, 2, 2)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/2   0]
            [  0   0 1/2]
            generated by 2 elements
            sage: OL = L.orthogonal_group()
            sage: f = OL.an_element()
            sage: fd = D.hom([d*f for d in D.smith_form_gens()])
            sage: OD(fd)
            [0 1 0]
            [1 0 0]
            [0 0 1]
        """
        from sage.groups.fqf_orthogonal.group import FqfOrthogonalGroup
        from sage.groups.fqf_orthogonal.group import _compute_gens
        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        from sage.groups.matrix_gps.linear import GL
        from sage.matrix.matrix_space import MatrixSpace
        from sage.rings.all import GF
        from copy import copy
        if not self.is_degenerate():
            gens = _compute_gens(self.normal_form()) #are we assuming normal form??
        elif self.invariants()[-1].is_prime():
            p = self.invariants()[-1]
            n = len(self.invariants())
            T = _normalize(self)
            q = T.gram_matrix_quadratic()
            k = len([i for i in range(n) if q[:,i]==0])
            r = n - k
            Idk = matrix.identity(k)
            Idr = matrix.identity(r)
            TR = T.submodule_with_gens(T.gens()[k:])
            gensTR = [TR._to_smith()*g*TR._to_gens() for g in _compute_gens(TR)]
            if k > 0:
                gens = [matrix.block_diagonal([Idk,g]) for g in gensTR]
                gens += [matrix.block_diagonal([g.matrix(),Idr]) for g in GL(k,p).gens()]
            else:
                gens = gensTR
            if k>0 and r>0:
                h = matrix.identity(n)
                for g in MatrixSpace(GF(p),r,k).basis():
                    h[k:,:k] = g
                    gens.append(copy(h))
            gens = [T._to_gens() * g * T._to_smith() for g in gens]
        else:
            raise NotImplementedError()
        ambient = AbelianGroupGap(self.invariants()).aut()
        gens = [ambient(g) for g in gens]
        gens = tuple(g for g in gens if g != ambient.one())
        return FqfOrthogonalGroup(ambient, gens, self)

    def primary_part(self, m):
        r"""
        Return the ``m``-primary part of this torsion quadratic module
        as a submodule.

        INPUT:

        - ``m`` -- an integer

        OUTPUT:

        - a submodule

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule((1/6)*ZZ^3,ZZ^3)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 6)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/36    0    0]
            [   0 1/36    0]
            [   0    0 1/36]
            sage: T.primary_part(2)
            Finite quadratic module over Integer Ring with invariants (2, 2, 2)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/4   0   0]
            [  0 1/4   0]
            [  0   0 1/4]

        TESTS::

            sage: T == T.primary_part(T.annihilator().gen())
            True
        """
        annihilator = self.annihilator().gen()
        a = annihilator.prime_to_m_part(m)
        return self.submodule( (a*self.V()).gens() )

    def submodule_with_gens(self, gens):
        r"""
        Return a submodule with generators given by ``gens``.

        INPUT:

        - ``gens`` -- a list of generators that convert into ``self``

        OUTPUT:

        - a submodule with the specified generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*10)
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: g = T.gens()
            sage: new_gens = [2*g[0], 5*g[0]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [2/5   0]
            [  0 1/2]

        The generators do not need to be independent::

            sage: new_gens = [g[0], 2*g[1], g[0], g[1]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10, 10)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/10    0 1/10    0]
            [   0  2/5    0  1/5]
            [1/10    0 1/10    0]
            [   0  1/5    0 1/10]

        TESTS:

        Test that things work without specified gens too::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
            sage: T.submodule(T.gens()[:2])
            Finite quadratic module over Integer Ring with invariants (5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0]
            [  0 1/5]
        """
        gens = [self(v) for v in gens]
        V = self.V().submodule([v.lift() for v in gens]) + self._W
        W = self.W()
        return TorsionQuadraticModule(V, W, gens=gens, modulus=self._modulus,
                                      modulus_qf=self._modulus_qf, check=False)

    def twist(self, s):
        r"""
        Return the torsion quadratic module with quadratic form scaled by ``s``.

        If the old form was defined modulo `n`, then the new form is defined
        modulo `n s`.

        INPUT:

        - ``s`` - a rational number

        EXAMPLES::

            sage: q = TorsionQuadraticForm(matrix.diagonal([3/9, 1/9]))
            sage: q.twist(-1)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 8/9]

        This form is defined modulo `3`::

            sage: q.twist(3)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/3Z:
            [  1   0]
            [  0 1/3]

        The next form is defined modulo `4`::

            sage: q.twist(4)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/4Z:
            [4/3   0]
            [  0 4/9]
        """
        s = self.base_ring().fraction_field()(s)
        n = self.V().degree()
        inner_product_matrix = s * self.V().inner_product_matrix()
        ambient = FreeQuadraticModule(self.base_ring(), n, inner_product_matrix)
        V = ambient.span(self.V().basis())
        W = ambient.span(self.W().basis())
        modulus = s * self._modulus
        modulus_qf = s * self._modulus_qf
        return TorsionQuadraticModule(V, W, modulus=modulus, modulus_qf=modulus_qf)

    def value_module(self):
        r"""
        Return `\QQ / m\ZZ` with `m = (V, W)`.

        This is where the inner product takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module()
            Q/Z
        """
        return QmodnZ(self._modulus)

    def value_module_qf(self):
        r"""
        Return `\QQ / n\ZZ` with `n\ZZ = (V,W) + \ZZ \{ (w,w) | w \in W \}`.

        This is where the torsion quadratic form takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module_qf()
            Q/2Z
        """
        return QmodnZ(self._modulus_qf)

    def to_smith(self):
        r"""
        Return the transformation matrix from the user to smith form generators.

        To go in the other direction use :meth:`to_gens`.

        OUTPUT:

        - an integer matrix

        EXAMLES::

            sage: L2 = IntegralLattice(3 * matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: D = L2.discriminant_group().normal_form()
            sage: D
            Finite quadratic module over Integer Ring with invariants (3, 6, 12)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0   0   0   0]
            [  0 1/4   0   0   0]
            [  0   0 1/3   0   0]
            [  0   0   0 1/3   0]
            [  0   0   0   0 2/3]
            sage: D.to_smith()
            [0 3 0]
            [0 0 3]
            [0 2 0]
            [1 0 0]
            [0 0 4]
            sage: T = D.to_smith()*D.to_gens()
            sage: T
            [ 3  0 33  0  0]
            [ 0 33  0  0  3]
            [ 2  0 22  0  0]
            [ 0  0  0  1  0]
            [ 0 44  0  0  4]

        The matrix `T` now satisfies a certain congruence::

            sage: for i in range(T.nrows()):
            ....:     T[:,i] = T[:,i] % D.gens()[i].order()
            sage: T
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        to_smith = matrix(self.base_ring(), [t.vector() for t in self.gens()])
        return to_smith

    @cached_method
    def to_gens(self):
        r"""
        Return the transformation matrix from smith form to user generators.

        To go in the other direction use :meth:`to_smith`.

        OUTPUT:

        - an integer matrix

        EXAMPLES::

            sage: L2 = IntegralLattice(3 * matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: D = L2.discriminant_group().normal_form()
            sage: D
            Finite quadratic module over Integer Ring with invariants (3, 6, 12)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0   0   0   0]
            [  0 1/4   0   0   0]
            [  0   0 1/3   0   0]
            [  0   0   0 1/3   0]
            [  0   0   0   0 2/3]
            sage: D.to_gens()
            [ 0  0  0  1  0]
            [ 1  0 11  0  0]
            [ 0 11  0  0  1]
            sage: T = D.to_gens()*D.to_smith()
            sage: T
            [ 1  0  0]
            [ 0 25  0]
            [ 0  0 37]

        This matrix satisfies the congruence::

            sage: for i in range(T.ncols()):
            ....:     T[:, i] = T[:, i] % D.smith_form_gens()[i].order()
            sage: T
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        invs = self.invariants()
        R = IntegerModRing(invs[-1])
        E = matrix.identity(R, len(invs))
        # view self as a submodule of (ZZ/nZZ)^(len(invs))
        B = self.to_smith().change_ring(R)
        for k in range(B.ncols()):
            B[:, k] *= invs[-1] // invs[k]
            E[:, k] *= invs[-1] // invs[k]
        to_gens = B.solve_left(E)
        to_gens = to_gens.change_ring(ZZ)
        # double check that we did not mess up
        # tmp = to_gens * self.to_smith()
        # E = matrix.identity(R, len(invs))
        # for k in range(E.ncols()):
        #     assert 0 == (tmp[:, k]-E[:, k]) % invs[k]
        return to_gens

    def direct_sum(self, other):
        r"""
        """
        V, fVs, fVo = self.V().direct_sum(other.V(),return_embeddings=True)
        W = self.W().direct_sum(other.W())
        n = len(self.V().gens())
        T = TorsionQuadraticModule(V, W, modulus=self._modulus)
        fs = self.hom([fVs(g.lift()) for g in self.gens()], T)
        fo = other.hom([fVo(g.lift()) for g in other.gens()], T)
        return T, fs, fo

    def _subgroup_to_gap(self, S):
        r"""
        """
        A = self.orthogonal_group().domain()
        gensS = [A(self(s)).gap() for s in S.gens()]
        Sgap = A.gap().Subgroup(gensS)
        return Sgap

    def _subgroup_from_gap(self, S):
        r"""
        """
        A = self.orthogonal_group().domain()
        S = self.submodule([self.linear_combination_of_smith_form_gens(
                            A(g).exponents())
                            for g in S.GeneratorsOfGroup()])
        return S

    def stab(self, G, S):
        r"""
        """
        Oq = self.orthogonal_group()
        mu = libgap.function_factory("mu:=function(x,g) return(Image(g,x)); end;")
        stab = libgap.Stabilizer(G.gap(), Sgap, mu)
        return stab

    def subgroup_representatives(self, H, G, algorithm="hulpke", return_stabilizers=False,
                                 order=None, max_order=None, g=None):
        r"""
        Return representatives of the subgroups of `H` modulo the action of `G`.

        INPUT:

        - ``H`` -- a submodule of `self`

        - ``G`` -- a group of automorphisms

          * ``"hulpke"`` -- following an algorithm of A. Hulpke

          * ``"brute force""`` -- enumerates all subgroups first and takes orbits.

        OUTPUT:

        - a list of submodules

        EXAMPLES::

            sage: T = TorsionQuadraticForm(matrix.diagonal([2/3,2/9,4/9]))
            sage: G = T.orthogonal_group()
            sage: subs1 = T.subgroup_representatives(T, G, algorithm="hulpke")
            sage: subs2 = T.subgroup_representatives(T, G, algorithm="brute force")
            sage: len(subs1) == len(subs2)
            True
        """
        if max_order is None:
            max_order = self.cardinality()
        A = G.domain()
        Hgap = A.subgroup([A(self(h)) for h in H.gens()]).gap()
        from sage.libs.gap.libgap import libgap
        if algorithm=="brute force":
            mu = libgap.function_factory("mu:=function(x,g) return(Image(g,x)); end;")
            all_subgroups = Hgap.AllSubgroups()
            subgroup_reps = G.gap().ExternalOrbitsStabilizers(all_subgroups, mu)
            # take a representative in each orbit
            subgroup_reps = [dict([['repr', S.Representative()],['stab', S.Stabilizer()]]) for S in subgroup_reps]
        elif algorithm == "hulpke":
            from sage.env import SAGE_EXTCODE
            gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits/subgroup_orbits.g'
            libgap.Read(gapcode)
            subgroup_reps = libgap.function_factory("SubgroupRepresentatives")
            subgroup_reps = [dict(S) for S in subgroup_reps(Hgap, G, max_order)]
        elif algorithm =="elementary abelian":
            from sage.env import SAGE_EXTCODE
            gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits/subgroup_orbits.g'
            libgap.Read(gapcode)
            subgroup_reps = libgap.function_factory("SubgroupRepresentatives_elementary_equiv")
            subgroup_reps = [dict(S) for S in subgroup_reps(Hgap, G, order, g)]
        else:
            raise ValueError("not a valid algorithm")
        # convert back to sage
        if return_stabilizers:
            subgroup_reps = [[self._subgroup_from_gap(S['repr']),G.subgroup(S['stab'].GeneratorsOfGroup())] for S in subgroup_reps]
        else:
            subgroup_reps = [self._subgroup_from_gap(S['repr']) for S in subgroup_reps]
        #subgroup_reps = map(self._subgroup_from_gap, subgroup_reps)
        return subgroup_reps

    def all_primitive_modulo(self, H1, H2, G, algorithm='hulpke'):
        r"""
        Return all totally isotropic subgroups `S` of `H1 + H2` such that
        ``H1 & S == 1`` and ``H2 & S = 1`` modulo the subgroup
        G of the orthogonal group of self

        Input:

        - ``H1``, ``H2`` -- subgroups of self

        - ``G`` - a subgroup of the automorphism group of self

        - ``algorithm`` -- a string which is either

          * ``"hulpke"`` -- following an algorithm of A. Hulpke

          * ``"brute force""`` -- enumerates all subgroups first and takes orbits.

        EXAMPLES::

            sage: q1 = TorsionQuadraticForm(matrix.diagonal([2/3,2/27]))
            sage: q2 = TorsionQuadraticForm(matrix.diagonal([2/3,2/9]))
            sage: q, fs , fo = q1.direct_sum(q2)
            sage: q.all_primitive_modulo(3*fs.image(),fo.image())
        """
        primitive_extensions = []
        A = G.domain()
        H = self.submodule(H1.gens() + H2.gens())
        if len(H.invariants()) > 6:
            print("this might take a while. Invariants: %s" %(H.invariants(),))
        subgroup_reps = self.subgroup_representatives(H, G, algorithm=algorithm)
        # filter for primitive and isotropic
        for S in subgroup_reps:
            if S.gram_matrix_quadratic() == 0:
                if S.V() & H1.V() == self.W() and S.V() & H2.V() == self.W():
                    primitive_extensions.append(S)
        return primitive_extensions

    def all_primitive_prime(self, other, H1, H2, G1, G2):
        r"""
        Return all totally isotropic subgroups `S` of `H1 + H2` such that
        ``H1 & S == 1`` and ``H2 & S = 1`` modulo the subgroup
        G of the orthogonal group of self

        Input:

        - ``H1``, ``H2`` -- subgroups of self

        - ``G`` - a subgroup of the automorphism group of self

        - ``algorithm`` -- a string which is either

          * ``"hulpke"`` -- following an algorithm of A. Hulpke

          * ``"brute force""`` -- enumerates all subgroups first and takes orbits.

        EXAMPLES::

            sage: q1 = TorsionQuadraticForm(matrix.diagonal([2/3,2/27]))
            sage: q2 = TorsionQuadraticForm(matrix.diagonal([2/3,2/9]))
            sage: q, fs , fo = q1.direct_sum(q2)
            sage: q.all_primitive_modulo(3*fs.image(),fo.image())
        """
        D, i1, i2 = self.direct_sum(other)
        primitive_extensions = []
        if H1.cardinality()==1 or H2.cardinality() == 1:
            return [self.submodule([])]
        p1 = H1.invariants()[-1]
        p2 = H2.invariants()[-1]
        if p1 != p2:
            raise ValueError("invariants do not match")
        if not p1.is_prime():
            raise ValueError("not a prime number")
        m = min(H1.cardinality(),H2.cardinality())
        subs1 = self.subgroup_representatives(H1, G1, return_stabilizers=True, max_order=m)
        subs2 = other.subgroup_representatives(H2, G2, return_stabilizers=True, max_order=m)
        #subs1 = [[_normalize(s[0]),s[1]] for s in subs1]
        subs2 = [[_normalize(s[0]),s[1]] for s in subs2]

        subs1 = [[s[0], s[0].orthogonal_group().subgroup(s[1].gens())] for s in subs1]
        subs2 = [[s[0], s[0].orthogonal_group().subgroup(s[1].gens())] for s in subs2]

        for S1 in subs1:
            S1n = _normalize(S1[0].twist(-1))
            if S1[0].cardinality()==1:
                primitive_extensions.append(D.submodule([]))
                continue
            n = len(S1[0].gens())
            for S2 in subs2:
                q1 = S1n.gram_matrix_quadratic()
                q2 = S2[0].gram_matrix_quadratic()
                if q1 != q2:
                    continue
                # there is a glue map

                A1 = S1[0].orthogonal_group().domain()
                A2 = S2[0].orthogonal_group().domain()
                gens1 = [A1(g).gap() for g in S1n.gens()]
                gens2 = [A2(g).gap() for g in S2[0].gens()]

                phi = A1.gap().GroupHomomorphismByImages(A2.gap(), gens1, gens2)
                O2 = S2[0].orthogonal_group()

                stab1phi= [phi.InducedAutomorphism(g) for g in S1[1].gap().GeneratorsOfGroup()]
                stab1phi = O2.gap().Subgroup(stab1phi)
                reps = O2.gap().DoubleCosetRepsAndSizes(S2[1],stab1phi)

                for g in reps:
                    g = g[0]
                    g = O2(g)
                    gens = [i1(S1n.gen(k)) + i2(S2[0].gen(k)*g) for k in range(n)]
                    primitive_extensions.append(D.submodule(gens))
        return primitive_extensions


    def all_primitive_prime_equiv(self, other, H1, H2, G1, G2, h1, h2, g):
        r"""
        Return all totally isotropic subgroups `S` of `H1 + H2` such that
        ``H1 & S == 1`` and ``H2 & S = 1`` modulo the subgroup
        G of the orthogonal group of self

        Input:

        - ``H1``, ``H2`` -- subgroups of ``self``, ``other``
        - ``G1``, ``G2`` -- subgroups of ``self``, ``other``
        - ``h1``, ``h2`` -- elements in the center of ``G1``, ``G2``
        - ``g`` -- non negative integer; the glue order `p^g`

        EXAMPLES::

            sage: q1 = TorsionQuadraticForm(matrix.diagonal([2/3,2/27]))
            sage: q2 = TorsionQuadraticForm(matrix.diagonal([2/3,2/9]))
            sage: q, fs , fo = q1.direct_sum(q2)
            sage: q.all_primitive_modulo(3*fs.image(),fo.image())
        """
        h1 = G1.subgroup([h1]).gen(0)
        h2 = G2.subgroup([h2]).gen(0)
        if not h1 in G1:#.center():
            raise ValueError()
        if not h2 in G2:#.center():
            raise ValueError()
        g = ZZ(g)
        if not g >= 0:
            raise ValueError()
        D, i1, i2 = self.direct_sum(other)
        OD = D.orthogonal_group()
        embedG1, embedG2 = direct_sum_embed(D, i1, i2, OD, G1, G2)

        from sage.libs.gap.libgap import libgap
        from sage.env import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits/subgroup_orbits.g'
        libgap.Read(gapcode)
        OnSubgroups = libgap.function_factory("OnSubgroups")

        primitive_extensions = []
        if g == 0:
            gens = embedG1.Image(G1.gap()).GeneratorsOfGroup()
            gens = gens.Concatenation(embedG2.Image(G2.gap()).GeneratorsOfGroup())
            return [[D.submodule([]),OD.subgroup(gens)]]
        if H1.cardinality()==1 or H2.cardinality()==1:
            return []
        p1 = H1.invariants()[-1]
        p2 = H2.invariants()[-1]
        if p1 != p2:
            raise ValueError("invariants do not match")
        if not p1.is_prime():
            raise ValueError("not a prime number")
        p = p1

        glue_order = p**g

        # these may not be invariant subspaces!!! ---> crap
        subs1 = self.subgroup_representatives(H1, G1, algorithm="elementary abelian",
                                              return_stabilizers=True, order=glue_order, g=h1.gap())
        subs2 = other.subgroup_representatives(H2, G2, algorithm="elementary abelian",
                                               return_stabilizers=True, order=glue_order, g=h2.gap())
        subs1 = [[_normalize(s[0].twist(-1)),s[0],s[1]] for s in subs1]
        subs2 = [[_normalize(s[0]),s[0],s[1]] for s in subs2]


        #subs1 = [[s[0], s[0].orthogonal_group().subgroup(s[1].gens())] for s in subs1]
        #subs2 = [[s[0], s[0].orthogonal_group().subgroup(s[1].gens())] for s in subs2]



        for S1 in subs1:
            S1n, S1, stab1 = S1
            n = len(S1.gens())
            for S2 in subs2:
                S2n, S2, stab2 = S2
                q1 = S1n.gram_matrix_quadratic()
                q2 = S2n.gram_matrix_quadratic()
                if q1 != q2:
                    continue
                # there is a glue map

                O1 = S1.orthogonal_group()
                O2 = S2.orthogonal_group()

                A1 = O1.domain()
                A2 = O2.domain()
                gens1 = [A1(g).gap() for g in S1n.gens()]
                gens2 = [A2(g).gap() for g in S2n.gens()]

                phi = A1.gap().GroupHomomorphismByImages(A2.gap(), gens1, gens2)

                h1_on_S2 = phi.InducedAutomorphism(O1(h1).gap())
                h2_on_S2 = O2(h2).gap()
                if not O2.gap().IsConjugate(h1_on_S2, h2_on_S2):
                    # this glue map cannot be modified to be equivariant
                    continue
                else:
                    # make it equivariant
                    g = O2.gap().RepresentativeAction(h2_on_S2, h1_on_S2)
                    phi = phi*g
                assert h2_on_S2==phi.InducedAutomorphism(O1(h1).gap())

                center = O2.gap().Centraliser(h2_on_S2)

                stab1phi= [phi.InducedAutomorphism(O1(g).gap()) for g in stab1.gens()]
                stab1phi = center.Subgroup(stab1phi)
                stab2c = O2.subgroup(stab2.gens())
                reps = center.DoubleCosetRepsAndSizes(stab2c,stab1phi)

                if reps.Size() > 1:
                    print("warning, more than one glue")
                for g in reps:
                    g = g[0]
                    phig = phi*g
                    g = O2(g)
                    gens = [i1(S1n.gen(k)) + i2(S2n.gen(k)*g) for k in range(n)]
                    ext = D.submodule(gens)
                    # we also need the centralizer of h1 x h2 in S1 x S2
                    ###############################
                    phig_graph = OD.domain().subgroup([OD.domain()(a) for a in gens]).gap()
                    S1_times_S2 =  [embedG1.Image(s.gap()) for s in stab1.gens()]
                    S1_times_S2 += [embedG2.Image(s.gap()) for s in stab2.gens()]
                    S1_times_S2 = OD.gap().Subgroup(S1_times_S2)
                    stab = S1_times_S2.Stabilizer(phig_graph, OnSubgroups).GeneratorsOfGroup()
                    ##############################
                    primitive_extensions.append([ext, OD.subgroup(stab)])
        return primitive_extensions




def _brown_indecomposable(q, p):
    r"""
    Return the Brown invariant of the indecomposable form ``q``.

    The values are taken from Table 2.1 in [Shim2016]_.

    INPUT:

    - ``q`` - an indecomposable quadratic form represented by a
      rational `1 \times 1` or `2 \times 2` matrix
    - ``p`` - a prime number

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import _brown_indecomposable
        sage: q = Matrix(QQ, [1/3])
        sage: _brown_indecomposable(q,3)
        6
        sage: q = Matrix(QQ, [2/3])
        sage: _brown_indecomposable(q,3)
        2
        sage: q = Matrix(QQ, [5/4])
        sage: _brown_indecomposable(q,2)
        5
        sage: q = Matrix(QQ, [7/4])
        sage: _brown_indecomposable(q,2)
        7
        sage: q = Matrix(QQ, 2, [0,1,1,0])/2
        sage: _brown_indecomposable(q,2)
        0
        sage: q = Matrix(QQ, 2, [2,1,1,2])/2
        sage: _brown_indecomposable(q,2)
        4
    """
    v = q.denominator().valuation(p)
    if p == 2:
        # brown(U) = 0
        if q.ncols() == 2:
            if q[0,0].valuation(2)>v+1 and q[1,1].valuation(2)>v+1:
                # type U
                return mod(0, 8)
            else:
                # type V
                return mod(4*v, 8)
        u = q[0,0].numerator()
        return mod(u + v*(u**2 - 1)/2, 8)
    if p % 4 == 1:
        e = -1
    if p % 4 == 3:
        e = 1
    if v % 2 == 1:
        u = q[0,0].numerator()//2
        if legendre_symbol(u,p) == 1:
            return mod(1 + e, 8)
        else:
            return mod(-3 + e, 8)
    return mod(0, 8)

def _normalize(D):
    r"""
    INPUT:

    - a possibly degenerate torsion quadratic form over a field

    OUTPUT:

    - a normal form
    """
    if D.cardinality() == 1:
        return D
    D = D.normal_form()
    if not D.is_degenerate():
        return D
    p = D.invariants()[-1]
    if not p.is_prime():
        raise ValueError("")
    gens = list(D.gens())
    kerb = D.orthogonal_submodule_to(D)
    n = len(gens)

    kergens = [gens[i] for i in range(n) if gens[i] in kerb]
    nondeg = [gens[i] for i in range(n) if not gens[i] in kerb]
    k = len(kergens)

    for i in range(k):
        if kergens[i].q()==1:
            for j in range(i+1, k):
                if kergens[j].q() == 1:
                    kergens[j] += kergens[i]
            kergens = kergens[:i] + kergens[i+1:] + [kergens[i]]
            break
    gens = kergens + nondeg
    # translate all others to the kernel of q
    return D.submodule_with_gens(gens)


def direct_sum_embed(D, i1, i2, OD, G1, G2):
    r"""

    EXAMPLES::

    sage: from sage.modules.torsion_quadratic_module import direct_sum_embed
    sage: T = TorsionQuadraticForm(matrix.diagonal([2/3]*3))
    sage: T1 = TorsionQuadraticForm(matrix.diagonal([2/3]*3))
    sage: T2 = TorsionQuadraticForm(matrix.diagonal([2/3]*2))
    sage: D, i1, i2 = T1.direct_sum(T2)
    sage: OD = D.orthogonal_group()
    sage: G1 = T1.orthogonal_group()
    sage: G2 = T2.orthogonal_group()
    sage: phi1, phi2 = direct_sum_embed(D,i1,i2,OD,G1,G2)
    """
    Dgap = OD.domain()
    direct_sum_gens = [(Dgap(D(g))).gap() for g in i1.image().gens() + i2.image().gens()]

    n1 = len(i1.image().gens())
    n2 = len(i2.image().gens())
    gensG1_imgs = []
    for f in G1.gens():
        imgs = [Dgap(i1(a*f)).gap() for a in i1.domain().gens()]
        imgs += direct_sum_gens[n1:]
        f = Dgap.gap().GroupHomomorphismByImages(Dgap.gap(), direct_sum_gens, imgs)
        gensG1_imgs.append(f)
    embed1 = G1.gap().GroupHomomorphismByImages(OD.gap(),[g.gap() for g in G1.gens()], gensG1_imgs)

    gensG2_imgs = []
    for f in G2.gens():
        imgs = [Dgap(i2(a*f)).gap() for a in i2.domain().gens()]
        imgs = direct_sum_gens[:n1] + imgs
        f = Dgap.gap().GroupHomomorphismByImages(Dgap.gap(), direct_sum_gens, imgs)
        gensG2_imgs.append(f)
    embed2 = G2.gap().GroupHomomorphismByImages(OD.gap(),[g.gap() for g in G2.gens()], gensG2_imgs)
    return embed1, embed2
