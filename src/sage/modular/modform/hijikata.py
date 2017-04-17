#########################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

r"""
Hijikata trace formula

Compute the trace `tr(T_n)` of the `n`th Hecke operator acting on
`S_k(\Gamma_0(N))`, for any `n \geq 1`, except if `n|N`, in which case
`n` must be prime.

AUTHORS:

- William A. Stein (February 5, 2012)
"""
from __future__ import print_function

from sage.rings.all import Integer, QQ, QuadraticField
from sage.functions.all import ceil, sign
from sage.misc.all import prod
from sage.libs.pari.all import pari
from sage.matrix.all import zero_matrix
from sage.arith.all import euler_phi, sigma
from sage.arith.misc import dedekind_psi
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.finite_rings.integer_mod_ring import Zmod


def w(d):
    r"""
    Return the `w` factor for the lattice of discriminant `d`.

    This is always 1 except when `d` is -3 or -4.

    See QuadraticField(d).number_of_roots_of_unit() / 2

    EXAMPLES::

        sage: from sage.modular.modform.hijikata import w
        sage: [w(d) for d in range(-7, -2)]
        [1, 1, 1, 2, 3]
    """
    if d == -4:
        return 2
    if d == -3:
        return 3
    return 1


def tof(a):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import tof
        sage: [tof(u) for u in [45, 49, 693]]
        [3, 7, 3]
    """
    a = Integer(a)
    t = prod(p ** ((e - (e % 2)) // 2) for p, e in a.factor())
    if Integer(a / t ** 2) % 4 == 1:
        return t
    return Integer(t / 2)


def sig(n, N):
    r"""
    This is a modified version of the sigma function that maps an integer n
    to the sum of the divisors of n.

    Values are modified only for divisors of N.

    See sage.rings.arith.sigma

    EXAMPLES::

        sage: from sage.modular.modform.hijikata import sig
        sage: [sig(d, 6) for d in range(1, 12)]
        [1, 2, 3, 7, 6, 6, 8, 15, 13, 18, 12]
    """
    N = Integer(N)
    n = Integer(n)
    if N % n == 0:
        return n
    return sigma(n, 1)
    # return prod((1 - p ** (e + 1)) / (1 - p) for p, e in n.factor())


def quadpoints(s, n, p, v):
    r"""
    Return the solutions of the quadratic equation `x ^ 2 - s x + n = 0`
    in the ring `\Zmod(p^v)`.

    .. TODO::

        find a better algorithm (using p-adics ?)

    EXAMPLES::

        sage: from sage.modular.modform.hijikata import quadpoints
        sage: quadpoints(8121,4471,691,1)
        [184, 336]
        sage: quadpoints(8121,4471,691,2)
        [219922, 265680]
    """
    p_pow = p
    sols = [x for x in range(p_pow) if (x ** 2 - s * x + n) % p_pow == 0]
    for k in range(v - 1):
        p_pow_next = p_pow * p
        sols_next = []
        for x in sols:
            for i in range(p):
                nx = x + i * p_pow
                if (nx ** 2 - s * nx + n) % p_pow_next == 0:
                    sols_next += [nx]
        p_pow = p_pow_next
        sols = sols_next
    return sols


def A(s, f, n, p, v):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import A
        sage: A(433,56,4,2833,1)
        2
    """
    rho = Integer(f).valuation(p)
    sr = set([x % p ** (v + rho) for x in quadpoints(s, n, p, v + 2 * rho)])
    return len([x for x in sr if
                (2 * x - s) % (p ** rho) == 0
                and ((n != p) or ((n == p) and x % p != 0))])


def B(s, f, n, p, v):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import B
        sage: B(433,56,4,191,1)
        2
    """
    rho = Integer(f).valuation(p)
    return len(set([x % p ** (v + rho)
                    for x in quadpoints(s, n, p, v + 2 * rho + 1)]))


def cp(s, f, n, p, v):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import cp
        sage: cp(434,2,4,79,1)
        2
    """
    ans = A(s, f, n, p, v)
    if Integer((s ** 2 - 4 * n) / f ** 2) % p == 0:
        ans += B(s, f, n, p, v)
    return ans


def c(s, f, n, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import c
        sage: c(434,2,4,79*23)
        4
    """
    return prod(cp(s, f, n, p, e) for p, e in N.factor())


def type_p(n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import type_p
        sage: type_p(471**2,4,191)
        104487111
    """
    if n.is_square():
        s = Integer(int((4 * n).sqrt()))
        return (Integer(1) / 4 * (s / 2) * n ** ((k // 2) - 1)
                * (c(s, 1, n, N) + (-1) ** k * c(-s, 1, n, N)))
    return 0


def absxy(s, n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import absxy
    """
    t = Integer(int((s ** 2 - 4 * n).sqrt()))
    x = (s - t) / 2
    y = (s + t) / 2
    product = (min(abs(x), abs(y)) ** (k - 1) / abs(x - y)) * sign(x) ** k
    return product * sum([euler_phi(Integer(t / f)) * c(s, f, n, N) / 2
                          for f in t.divisors()])


def type_h(n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import type_h
        sage: type_h(471**2,4,191)
        7741300
    """
    start = int(ceil(2 * n.sqrt()))
    if n.is_square():
        start += 1
    return sum([QQ(absxy(s, n, k, N) + absxy(-s, n, k, N)) for
                s in range(start, 4 * n + 1) if (s ** 2 - 4 * n).is_square()])


def classno(d, proof=True):
    r"""
    Return the class number of the order of discriminant d.

    See QuadraticField(d).class_number()

    EXAMPLES::

        sage: from sage.modular.modform.hijikata import classno
        sage: classno(-163)
        1
    """
    f = 1 if proof else 0
    return Integer(pari(d).qfbclassno(flag=f))


def xy(s, n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import xy
        sage: xy(4,59,6,7)
        7240
        sage: xy(81,14,3,7)
        162
    """
    K = QuadraticField(s ** 2 - 4 * n)
    a = K.gen()
    # x and y are the solutions to X^2 - s*X + n = 0.
    x = (s + a) / 2
    y = (s - a) / 2
    product = Integer(1) / 2 * (x ** (k - 1) - y ** (k - 1)) / (x - y)
    return QQ(product * sum(classno(Integer((s ** 2 - 4 * n) / f ** 2))
                            / w((s ** 2 - 4 * n) / f ** 2) * c(s, f, n, N)
                            for f in tof(s ** 2 - 4 * n).divisors()))


def type_e(n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import type_e
        sage: type_e(47**2,4,19)
        227310
    """
    r = int(2 * n.sqrt())
    if n.is_square():
        r -= 1
    # return sum([xy(s, n, k, N) for s in range(-r, r+1)])
    # WARNING: I *conjecture* that the sum below is the same as the one
    # that I commented out above.
    return xy(0, n, k, N) + 2 * sum(xy(s, n, k, N) for s in range(1, r + 1))


def type_e_conj(n, k, N):
    """
    Conjectural formula for type_e

    EXAMPLES::

        sage: from sage.modular.modform.hijikata import type_e_conj
        sage: type_e_conj(47**2,4,19)
        227310
    """
    r = int(2 * n.sqrt())
    if n.is_square():
        r -= 1
    return sum([xy(s, n, k, N) for s in range(-r, r+1)])


def sum_s(n, k, N):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import sum_s
        sage: sum_s(47**2,4,19)
        331135
    """
    return type_p(n, k, N) + type_h(n, k, N) + type_e(n, k, N)


def test_trace_hecke_operator(n_range, k_range, N_range, verbose=False):
    r"""
    Verify that the trace_hecke_operator command gives the same output
    as the trace computed using modular symbols.  If not, a RuntimeError
    exception is raised.

    INPUT:

    - ``n_range`` -- list of indexes `n` for which the Hecke operator `T_n`
      is computed when `n` is either prime or `n` is coprime to the level
    - ``k_range`` -- list of weights; odd weights are ignored (since
      trace is always 0)
    - ``N_range`` -- list levels
    - ``verbose`` -- bool (default: ``True``); if ``True`` print level and
      weight

    EXAMPLES:

    Test that level 1 traces are correct::

        sage: hijikata.test_trace_hecke_operator([1..14],[2, 4, .., 40],[1]) # long time

    Test levels up to 20 and weights 2, 4::

        sage: hijikata.test_trace_hecke_operator([1..14],[2,4],[1..20],) # long time
    """
    from sage.modular.modsym.all import ModularSymbols
    for N in N_range:
        for k in k_range:
            if k % 2 == 0:
                if verbose:
                    print("(N, k) =", (N, k), end="")
                S = ModularSymbols(N, k, sign=1).cuspidal_submodule()
                for n in n_range:
                    n = Integer(n)
                    if n.gcd(N) == 1 or n.is_prime():
                        if S.hecke_operator(n).trace() != trace_hecke_operator(n, k, N):
                            raise RuntimeError("trace_hecke_operator(n=%s,k=%s,N=%s) disagrees with modular symbols trace" % (n, k, N))
                if verbose:
                    print(" (good)")


def trace_hecke_operator(n, k, N=1):
    r"""
    Return the trace of the Hecke operator `T_n` on
    `S_k(\Gamma_0(N)))`.  The only constraint on `n` is that if
    `GCD(n,N) \neq 1`, then `n` must be prime.

    INPUT:

    - `n` -- positive integer
    - `k` -- integer at least 2
    - `N` -- positive integer

    OUTPUT:

    - rational number

    EXAMPLES:

    We compute the trace of the first few Hecke operators on level 1
    weight 12::

        sage: hijikata.trace_hecke_operator(1, 12)
        1
        sage: v = [hijikata.trace_hecke_operator(n, 12) for n in [1..30]]
        sage: v[:10]
        [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]
        sage: list(delta_qexp(31))[1:] == v
        True

    We can compute dimensions using the trace formula::

        sage: hijikata.trace_hecke_operator(1, 20000)
        1666
        sage: dimension_cusp_forms(1, 20000)
        1666
        sage: hijikata.trace_hecke_operator(1, 2000, 11)
        1998
        sage: dimension_cusp_forms(11, 2000)
        1998

    An advantage of the trace formula is that can quickly compute
    traces of operators that one could not do using other techniques
    like modular symbols.  For example, here we instantly compute the
    trace of `T_2` on `S_{20000}(\Gamma_0(19))`::

        sage: t = hijikata.trace_hecke_operator(2, 20000, 19)
        sage: t.valuation(2)
        1
    """
    n = Integer(n)
    k = Integer(k)
    N = Integer(N)
    if n.gcd(N) != 1 and not n.is_prime():
        raise ValueError("n must be prime when n and N are not coprime")

    if k % 2 != 0:
        return 0
    if k == 2:
        t = sig(n, N)
    else:
        t = 0
    t += -sum_s(n, k, N)
    if n.is_square():
        t += (k - 1) * dedekind_psi(N) / 12 * n ** ((k // 2) - 1)
    return t


def trace_modular_form(k, prec):
    r"""
    Return the trace modular form sum `Tr(T_n) q^n` to absolute
    precision prec, where `T_n` is the `n`th Hecke operator on
    `S_k(SL_2(Z))`.

    The complexity is almost entirely a function of ``prec``, but not of
    `k`.

    INPUT:

    - `k` -- positive integer
    - ``prec`` -- positive integer

    EXAMPLES::

        sage: hijikata.trace_modular_form(12, 6)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
        sage: hijikata.trace_modular_form(24, 6)
        2*q + 1080*q^2 + 339480*q^3 + 25326656*q^4 + 73069020*q^5 + O(q^6)
    """
    R = QQ[['q']]
    return R([0] + [trace_hecke_operator(n, k, 1) for n in range(1, prec)],
             prec)


def Tp(p, r, k, f):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import Tp
        sage: q = PowerSeriesRing(QQ,'q').gen()
        sage: hijikata.Tp(3, 2, 4, 1+q+q**3+q**7+O(q**266))
        -27 + 27*q^3 + 729*q^9 + 729*q^27 + O(q^30)
        sage: hijikata.Tp(1,2,4,1+q+q**3+q**7+O(q**266))
        -1 + 3*q + 3*q^3 + 3*q^7 + O(q^266)
    """
    if r > 1:
        return (Tp(p, 1, k, Tp(p, r - 1, k, f))
                - p ** (k - 1) * Tp(p, r - 2, k, f))
    if r == 1:
        R = QQ[['q']]
        q = R.gen()
        prec = int(((f.prec() - 1) / p) + 1)
        return R(sum(f[n * p] * q ** n + p ** (k - 1) * f[n] * q ** (n * p)
                     for n in range(1, prec)), prec)
    if r == 0:
        return f


def hecke_operator(n, k, f):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import hecke_operator
        sage: q = PowerSeriesRing(QQ,'q').gen()
        sage: hecke_operator(8, 12, q+q**2+O(q**100))
        4194304*q^4 + 8589934592*q^8 + O(q^13)
    """
    n = Integer(n)
    for p, e in n.factor():
        f = Tp(p, e, k, f)
    return f


def trace_formula_basis(k, prec):
    r"""
    EXAMPLES::

        sage: hijikata.trace_formula_basis(24, 6)
        [2*q + 1080*q^2 + 339480*q^3 + 25326656*q^4 + 73069020*q^5 + O(q^6),
        1080*q + 42103872*q^2 + O(q^3)]

    Double check the claimed result above to higher precision::

        sage: B = [vector(QQ,f.qexp(10)) for f in CuspForms(1, 24).basis()]
        sage: C = [vector(QQ, f.padded_list(10)) for f
        ....: in hijikata.trace_formula_basis(24, 21)]
        sage: span(B) == span(C)
        True
    """
    f = trace_modular_form(k, prec)
    d = Integer(f[1])  # the dimension
    return [hecke_operator(n, k, f) for n in range(1, d + 1)]


def basis_matrix(B):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import basis_matrix
        sage: basis_matrix([[4,5,6],[6,7,8]])
        [5 6]
        [7 8]
    """
    d = len(B)
    I = zero_matrix(QQ, d)
    for r in range(d):
        for c in range(d):
            I[r, c] = B[r][c + 1]
    return I


def hecke_operator_matrix(k, n):
    r"""
    EXAMPLES::

        sage: hijikata.hecke_operator_matrix(24, 2)
        [       0 20468736]
        [       1     1080]
        sage: hijikata.hecke_operator_matrix(24, 2).fcp()
        x^2 - 1080*x - 20468736
        sage: CuspForms(1, 24).hecke_polynomial(2)
        x^2 - 1080*x - 20468736
    """
    from sage.modular.dims import dimension_cusp_forms
    d = dimension_cusp_forms(1, k)
    B = trace_formula_basis(k, n * d ** 2 + 1)
    return hecke_operator_matrix_wrt_basis(k, n, B)


def hecke_operator_matrix_wrt_basis(k, n, B):
    r"""
    EXAMPLES::

        sage: from sage.modular.modform.hijikata import hecke_operator_matrix_wrt_basis
    """
    d = len(B)
    I = basis_matrix(B) ** (-1)
    A = []
    for j in range(d):
        g = hecke_operator(n, k, B[j])
        A.append([g[i] for i in range(1, d + 1)])
    A = I.parent()(A)
    return I * A
