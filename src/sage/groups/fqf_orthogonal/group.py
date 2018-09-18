r"""
Orthogonal groups of torsion quadratic forms.

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Simon Brandhorst (2018-05-15): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.abelian_gps.abelian_aut import AbelianGroupAutomorphismGroup_subgroup, AbelianGroupAutomorphism, AbelianGroupAutomorphismGroup_gap
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
from sage.groups.fqf_orthogonal.gens import _gens
from sage.misc.cachefunc import cached_method
from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
from sage.rings.all import mod, ZZ, IntegerModRing, Zp, QQ
from sage.matrix.all import matrix
from sage.categories.action import Action

class FqfIsometry(AbelianGroupAutomorphism):
    r"""
    Isometry of a finite quadratic/bilinear form.

    INPUT:

    - ``parent`` -- the parent :class:`~FqfOrthogonalGroup`
    - ``x`` -- a libgap element
    - ``check`` -- bool (default: ``True``)
    """

    def _repr_(self):
        r"""
        Return the string represenation of ``self``.

        EXAMPLES::

            sage: q = matrix.diagonal([1/3,2/3])
            sage: q = TorsionQuadraticForm(q)
            sage: G = q.orthogonal_group()
            sage: g = G.gens()[0]
            sage: g
            [2 0]
            [0 2]
        """
        return str(self.matrix())

    def __call__(self, x):
        r"""
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: q = matrix.diagonal([1/3,2/3])
            sage: q = TorsionQuadraticForm(q)
            sage: G = q.orthogonal_group()
            sage: g = G.gens()[0]
            sage: g
            [2 0]
            [0 2]
            sage: x = q.gens()[0]
            sage: g(x)
            (2, 0)
        """
        if x in self.parent().invariant_form():
            return x*self
        else:
            return AbelianGroupAutomorphism.__call__(self, x)

    def _gap_init(self):
        r"""
        """
        s = self.gap().__repr__()
        s = s.replace("->", ",")
        s = "GroupHomomorphismByImages(A, A, " + s + ")"
        return s

class FqfOrthogonalGroup(AbelianGroupAutomorphismGroup_subgroup):
    r"""
    Return the orthogonal group of this torsion quadratic form.

    Do not call this class directly instead use
    :meth:`sage.modules.torsion_quadratic_module.orthogonal_group`.

    INPUT:

        - ``T`` a non degenerate torsion quadratic module.

    EXAMPLES::

        sage: q = matrix.diagonal(QQ, [2/3, 2/3, 4/3])
        sage: T = TorsionQuadraticForm(q)
        sage: T
        Finite quadratic module over Integer Ring with invariants (3, 3, 3)
        Gram matrix of the quadratic form with values in Q/2Z:
        [4/3   0   0]
        [  0 2/3   0]
        [  0   0 2/3]
        sage: T.orthogonal_group()
        Group of isometries of
        Finite quadratic module over Integer Ring with invariants (3, 3, 3)
        Gram matrix of the quadratic form with values in Q/2Z:
        [4/3   0   0]
        [  0 2/3   0]
        [  0   0 2/3]
        generated by 2 elements
        sage: q = matrix(QQ, 3, 3, [3/2, 0, 0, 0, 1/4, 0, 0, 0 , 1/4])
        sage: T = TorsionQuadraticForm(q)
        sage: T
        Finite quadratic module over Integer Ring with invariants (2, 4, 4)
        Gram matrix of the quadratic form with values in Q/2Z:
        [3/2   0   0]
        [  0 1/4   0]
        [  0   0 1/4]
        sage: T.orthogonal_group().order()
        8
    """
    Element = FqfIsometry

    def __init__(self, ambient, gens, fqf, check=True):
        r"""
        """
        # We act on the smith form generators
        # because they are independent
        if not isinstance(fqf, TorsionQuadraticModule):
            raise TypeError("input must be a torsion quadratic module")
        if not isinstance(ambient, AbelianGroupAutomorphismGroup_gap):
            raise TypeError("input must be a torsion quadratic module")
        if not fqf.invariants() == ambient.domain().gens_orders():
            raise ValueError("invariants of the abelian group do not match")
        gens = [ambient(g) for g in gens]
        self._invariant_form = fqf
        AbelianGroupAutomorphismGroup_subgroup.__init__(self, ambient, gens)
        if check:
            for g in self.gens():
                if not self._preserves_form(g):
                    raise ValueError("%s does not preserve the quadratic form"%g)

    def invariant_form(self):
        r"""
        """
        return self._invariant_form

    def _gap_init(self):
        r"""
        """
        invs = self.domain().gap().AbelianInvariants().__repr__()
        s = "A:=AbelianGroup(" + invs + ");\n"
        pcgs = self.domain().gap().Pcgs()
        s += "f:=Pcgs(A);"
        for n in range(1,len(pcgs)+1):
            s += "f%s:=f[%s];"%(n,n)
        s += "\n"
        s += "aut:= AutomorphismGroup(A);\n"
        s += "G:=Subgroup(aut, ["
        for g in self.gens():
            s += g._gap_init() + ",\n"
        s = s[:-1]
        s += "]);"
        return s

    def _element_constructor_(self, x, check=True):
        r"""
        Construct an element from ``x`` and handle conversions.

        INPUT:

        - ``x`` -- something that converts in can be:

          * a libgap element
          * an integer matrix in the covering matrix ring
          * a class:`sage.modules.fg_pid.fgp_morphism.FGP_Morphism`
            defining an automorphism -- the domain of ``x`` must have
            invariants equal to ``self.domain().gens_orders()``
            if x in self.invariant_form().W().orthogonal_group()
          * something that acts on the invariant form module

        EXAMPLES::

            sage: L = IntegralLattice("A4").twist(2).direct_sum(IntegralLattice("A2"))
            sage: D = L.discriminant_group()
            sage: OL = L.orthogonal_group()
            sage: OD = D.orthogonal_group()
            sage: f = OL.an_element()
            sage: fbar = OD(f)
            sage: fbar
            [ 0  0  0 15]
            [ 1  1  0  0]
            [ 1  0  0  0]
            [ 1  0  1  1]

        TESTS::

            sage: all([x*f==x*fbar for x in D.gens()])
            True
        """
        # the super class knows what to do
        from sage.libs.gap.element import GapElement
        if not type(x) is GapElement:
            try:
                # if there is an action try that
                gen = self.invariant_form().smith_form_gens()
                x = matrix(ZZ, [(g*x).vector() for g in gen])
            except TypeError:
                pass
        f = AbelianGroupAutomorphismGroup_subgroup._element_constructor_(self, x, check=check)
        if check:
            # check that the form is preserved
            # this is expensive
            if not self._preserves_form(f):
                raise ValueError("not an isometry")
        return f

    def _preserves_form(self, f):
        r"""
        Return if f preserves the form.
        """
        g = self.invariant_form().smith_form_gens()
        for i in range(len(g)):
            if (g[i]*f).q() != g[i].q():
                return False
            for j in range(i+1, len(g)):
                if (g[i]*f).b(g[j]*f) != (g[i]*f).b(g[j]*f):
                    return False
        return True

    def _get_action_(self, S, op, self_on_left):
        r"""
        Provide the coercion system with an action.

        EXAMPLES::

            sage: q = matrix.diagonal([1/3,2/3])
            sage: q = TorsionQuadraticForm(q)
            sage: G = q.orthogonal_group()
            sage: G._get_action_(q, operator.mul, False)
            Right action by Group of isometries of
            Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 1/3]
            generated by 2 elements on Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 1/3]
        """
        import operator
        if op == operator.mul and not self_on_left:
            T = self.invariant_form()
            if S == T:
                return ActionOnFqf(self, S)
            try:
                if S.is_submodule(T):
                    # check if the submodule is invariant
                    if all([T(s)*g in S for s in S.gens() for g in self.gens()]):
                        return ActionOnFqf(self, S, on_submodule=True)
            except AttributeError:
                pass
            try:
                return AbelianGroupAutomorphismGroup_subgroup._get_action_(self, S, op, self_on_left)
            except AttributeError:
                pass

    def _subgroup_constructor(self, libgap_subgroup):
        r"""
        Create a subgroup from the input.

        See :class:`~sage.groups.libgap_wrapper`. Override this in derived
        classes.

        EXAMPLES::

            sage: q = TorsionQuadraticForm(matrix.diagonal([2/3,2/3]))
            sage: G = q.orthogonal_group()
            sage: G.subgroup(G.gens()[:1])
            Group of isometries of
            Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/2Z:
            [2/3   0]
            [  0 2/3]
            generated by 1 elements
        """
        generators = libgap_subgroup.GeneratorsOfGroup()
        generators = tuple(self(g) for g in generators)
        return FqfOrthogonalGroup(self, generators, self.invariant_form())

    def _repr_(self):
        r"""
        The string represenation of ``self``.

        EXAMPLES::

            sage: q = TorsionQuadraticForm(matrix.diagonal([2/3,2/3]))
            sage: q.orthogonal_group()
            Group of isometries of
            Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/2Z:
            [2/3   0]
            [  0 2/3]
            generated by 2 elements
        """
        return "Group of isometries of \n%s\ngenerated by %s elements"%(self.invariant_form(), len(self.gens()))

    def det_spin_homomorphism(self, rank, det):
        r"""
        Return the det spin homomorphism.

        This has only sense if ``self`` is the orthogonal group of a
        discriminant group of a lattice.

        EXAMPLES::

            sage: L = IntegralLattice("A2").twist(5)
            sage: L = L.direct_sum(L.twist(-1))
            sage: Oq = L.discriminant_group().orthogonal_group()
            sage: Oq.det_spin_homomorphism().image(Oq).order()
            4
        """
        from sage.groups.fqf_orthogonal.spinor import GammaA
        S = (2*self.domain().order()).prime_divisors()
        T = self.invariant_form()
        L = T.W()
        Gamma = GammaA(S)
        sigma_sharp = Gamma.sigma_sharp(rank, det, T)
        codom = Gamma.quotient(sigma_sharp)

        from sage.groups.fqf_orthogonal.lift import Hensel_qf
        from sage.rings.all import Zp, QQ
        from sage.quadratic_forms.all import QuadraticForm
        det_spin = dict([[f,[]] for f in self.gens()])

        for p in S:
            Tp = T.primary_part(p).normal_form()
            if Tp.cardinality() == 1:
                for f in self.gens():
                    det_spin[f].append((p, (1, 1)))
                continue
            u = matrix([b.vector() for b in Tp.gens()])
            Op = Tp.orthogonal_group()
            q = Tp.gram_matrix_quadratic()
            # Shimada 5.2 Step 2 page 25
            # The discriminant group has full length and a term 1/2 or 3/2
            # to determine the lattice one needs to know the determinant
            val, unit = (q.det()*det).val_unit(ZZ(p))
            assert val==0
            if p == 2 and len(Tp.invariants())== rank and unit % 8 != 1:
                n = q.ncols()
                for i in range(1,n-1):
                    if q[i,i+1]==0 and q[i-1,i]==0 and q[i,i].denominator()==2:
                        q[i,i] *= 5
                        break
                else:
                    if q[0,0].denominator()==2 and (n==1 or q[0,1]==0):
                        q[0,0] *= 5
                    elif q[-1,-1].denominator()==2 and q[-1,-2]==0:
                        q[-1,-1] *= 5
                    else:
                        raise AssertionError('bug in det_spin_homomorphism')
            M = q.inverse()
            # diagonalize
            qf = QuadraticForm(QQ, M)
            diag, t = qf.rational_diagonal_form(return_matrix=True)
            diag = diag.Hessian_matrix()
            t = t.T
            q0 = q*q.denominator()
            v1 = t.denominator().valuation(p)
            v2 = t.inverse().denominator().valuation(p)
            v = -v1 -v2 # lower bound for precision loss due to diagonalization
            # compute the spin of lifts of the generators
            for f in self.gens():
                # take only the action on the p-part
                fp = Op(f).matrix()
                prec0 = 1
                prec = 25 # initial precision
                # change to the user basis
                g = u * fp * u.inverse()
                while True:
                    R = Zp(p, type='fixed-mod', prec=prec+3, print_mode='terse',
                        show_prec=False, print_pos=False)
                    g = Hensel_qf(q0.change_ring(R), g.change_ring(R), prec0, prec)
                    g = g.change_ring(ZZ)
                    try:
                        gg = t*M*g*M.inverse()*t.inverse()
                        det_p, spin_p = det_spin_p(diag, gg, p, prec + v)
                        det_spin[f].append((p, (det_p, spin_p)))
                        break
                    except ArithmeticError:
                        # retry with higher precision
                        prec0 = prec
                        prec = 2 * prec
        imgs = [Gamma.prod([Gamma.embed(det=ds[1][0], spin=ds[1][1], p=ds[0])
                            for ds in det_spin[f]])
                for f in self.gens()]
        return self.hom([codom(g) for g in imgs],codom)

def _compute_gens(T):
    r"""
    Return generators of the orthogonal group of ``T``.

    INPUT:

    - ``T`` -- torsion orthogonal module in normal form.

    OUTPUT:

    - a list of matrices -- the generators

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.group import _compute_gens
        sage: T = TorsionQuadraticForm(matrix.diagonal([2/3,2/3]))
        sage: _compute_gens(T.normal_form())
        [
        [  1 726]  [0 1]  [0 1]
        [  3   1], [2 0], [1 0]
        ]
    """
    assert T==T.normal_form()
    # corner case
    invs = T.invariants()
    if len(invs) == 0:
        return []

    # normal form gens for the different primes
    blocks = []
    gensT_orders = [t.order() for t in T.gens()]
    n = len(T.gens())
    P = T.annihilator().gen().prime_divisors()
    for p in P:
        indices = []
        for k in range(len(gensT_orders)):
            if mod(gensT_orders[k], p) == 0:
                indices.append(k)
        blocks.append([p, tuple(indices)])

    to_smith = T.to_smith()
    to_normal = T.to_gens()

    # compute generators of the orthogonal groups
    gens = []
    for bl in blocks:
        # compute the generators of the p-primary part
        # the whole group is the direct product of the p-primary parts
        p = bl[0]
        indices = bl[1]
        q_p = T.gram_matrix_quadratic()[indices, indices]
        b = invs[-1].valuation(p)
        G_p = q_p * p**b
        if p != 2:
            # make sure each homogeneous block of G_p stays in normal form
            G_p = G_p/2

        R = Zp(p, type='fixed-mod', prec=b+5, show_prec=False)
        G_p = G_p.change_ring(R)

        # the generators in matrix form
        gens_mat = _gens(G_p, b)

        # extend as identity on the orthogonal complement
        E1 = matrix.identity(indices[0])
        E2 = matrix.identity(n - indices[-1] - 1)
        for g in gens_mat:
            g = g.change_ring(ZZ)
            g = matrix.block_diagonal([E1,g,E2])
            g = to_normal * g * to_smith
            gens.append(g)
    return gens

class ActionOnFqf(Action):
    r"""
    Action on a finite quadratic module.

    INPUT:

    - ``orthogonal_grp`` --  an instance of :class:`GroupOfIsometries`
    - ``fqf`` -- a torsion quadratic module
    - ``on_submodule`` -- bool (default: ``False``)
    - ``is_left`` -- bool (default: ``False``)

    EXAMPLES::

        sage: q = matrix.diagonal([1/3,2/3])
        sage: q = TorsionQuadraticForm(q)
        sage: G = q.orthogonal_group()
        sage: g = G.gens()[0]
        sage: g
        [2 0]
        [0 2]
        sage: x = q.gens()[0]
        sage: x*g
        (2, 0)
    """
    def __init__(self, orthogonal_grp, fqf, on_submodule=False, is_left=False):
        r"""
        Initialize the action

        TESTS::

            sage: from sage.groups.fqf_orthogonal.group import ActionOnFqf
            sage: q = matrix.diagonal([1/3, 2/3])
            sage: q = TorsionQuadraticForm(q)
            sage: G = q.orthogonal_group()
            sage: A = ActionOnFqf(G, q, is_left=True)
            Traceback (most recent call last):
            ...
            ValueError: the action is from the right
        """
        import operator
        self._on_submodule = on_submodule
        if is_left:
            raise ValueError("the action is from the right")
        Action.__init__(self, orthogonal_grp, fqf, is_left, operator.mul)

    def _call_(self, a, g):
        r"""
        This defines the group action.

        INPUT:

        - ``a`` -- an element of the invariant submodule
        - ``g`` -- an element of the acting group

        OUTPUT:

        - an element of the invariant submodule

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.group import ActionOnFqf
            sage: q = matrix.diagonal([1/3, 2/3])
            sage: q = TorsionQuadraticForm(q)
            sage: G = q.orthogonal_group()
            sage: g = G.gens()[0]
            sage: A = ActionOnFqf(G, q)
            sage: A
            Right action by Group of isometries of
            Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 1/3]
            generated by 2 elements on Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 1/3]
            sage: x = q.an_element()
            sage: g = G.an_element()
            sage: A(x,g).parent()
            Finite quadratic module over Integer Ring with invariants (3, 3)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 1/3]
            sage: q = TorsionQuadraticForm(matrix.diagonal([1/3,1/3,3/8,1/4]))
            sage: G = q.orthogonal_group()
            sage: q2 = q.primary_part(2)
            sage: q2.gen(0)*G.gen(0)
            (3, 4)
        """
        if self.is_left():
            pass
            # this would be a left action but... we do not allow it.
            # v = (a.vector()*g.matrix().inverse())
            # P = a.parent()
            # return P.linear_combination_of_smith_form_gens(v)
        elif self._on_submodule:
            S = a.parent()
            T = g.parent().invariant_form()
            return S(T(a)*g)
        else:
            v = (a.vector()*g.matrix())
            P = a.parent()
            return P.linear_combination_of_smith_form_gens(v)


def det_spin_p(G, T, p, nu):
    r"""
    Return approximations for  `(det_p, spin_p)` of ``T``.

    The algorithm is due to Shimada.
    We follow the conventions of Miranda and Morrison that the quadratic form is defined by
    Q(x) = (x G x.T)/2. Then the spinor norm of the reflection in x is Q(x).

    INPUT:

        - ``G`` -- a diagonal matrix
        - ``T`` -- an isometry up to some precision
        - ``p`` -- a prime number
        - ``nu`` -- an integer giving the valuation of the approximation
        error of ``T``

    EXAMPLES::
    """
    from sage.groups.fqf_orthogonal.lift import _min_val
    from sage.rings.all import Qp
    def mv(A):
        return _min_val(A.change_ring(Qp(p)))
    if p == 2:
        delta = 1
    else:
        delta = 0
    gammaL = [d.valuation(p) for d in G.diagonal()]
    gamma = min(gammaL)
    l = G.ncols()
    E = G.parent()(1)
    reflection_vectors = []

    k = 0
    while k < l:
        g = T.row(k)
        # error estimates
        lambd = mv(g)
        rho = min(delta + nu + gamma, 2*nu + gamma)
        sigma = min( delta + nu + gamma, delta + nu + lambd, 2*nu + gamma)
        kappa = sigma - gammaL[k] - 2*delta
        if (rho <= gammaL[k] + delta) or (kappa < 1 + 2*delta):
            raise ArithmeticError("Recompute with higher precision") #or a ValueError ?
        bm = g - E.row(k)
        qm = bm * G * bm
        if qm.valuation(p) <= gammaL[k] + 2*delta:
            tau1 = reflection(G, bm)
            reflection_vectors.append(bm)
            tau2 = E
        else:
            bp = g + E.row(k)
            qp = bp * G * bp
            assert qp.valuation(p) <= gammaL[k] + 2*delta
            tau1 = reflection(G, bp)
            tau2 = reflection(G, E.row(k))
            reflection_vectors.append(bp)
            reflection_vectors.append(E.row(k))
        lambdaT = mv(T)
        alpha = mv(tau1)
        beta = mv(tau2)
        theta = gamma + min(kappa + 2*min(0,lambd), nu + min(0,lambd), 2*nu)
        nu = min(nu + alpha, lambdaT + theta - gammaL[k] - delta, nu + theta - gammaL[k] - delta) + beta
        T = T * tau1 * tau2
        k += 1
    err = mv(T-E)
    assert err >= nu
    spinor_norm = QQ.prod([v*G*v/2 for v in reflection_vectors])
    determinant = QQ(-1)**(len(reflection_vectors))
    v, u = spinor_norm.val_unit(p)
    if p == 2:
        u = u % 8
    else:
        u = u % p
    spinor_norm = u * p**(v % 2)
    return determinant, spinor_norm

def reflection(G, v):
    r"""
    Return the matrix represenation of the orthogonal reflection in `v`.

    INPUT:

    - ``v`` -- a vector
    - ``G`` -- a symmetric matrix

    OUTPUT:

    - a matrix

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.group import reflection
        sage: G = matrix.identity(2)
        sage: v = vector([1,0])
        sage: reflection(G, v)
        [-1  0]
        [ 0  1]
    """
    n = v.degree()
    E = v.parent().basis()
    vsq = v*G*v
    ref = []
    for k in range(n):
        tk = E[k] - 2/vsq*(E[k]*G*v)*v
        ref.append(tk)
    return matrix(ref)
