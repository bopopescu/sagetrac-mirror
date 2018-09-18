"Genus"

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#                          Gabriele Nebe <nebe@math.rwth-aachen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.misc.all import prod
from sage.arith.all import LCM
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.rational_field import RationalField, QQ
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from copy import copy, deepcopy
from sage.misc.misc import verbose

def all_genera_by_det(sig_pair, determinant, max_scale=None, even=True):
    r"""
    Return a list of all global genera with the given conditions.

    Here a genus is called global if it is non empty.

    INPUT:

    - ``sig_pair`` -- a pair of non-negative integers giving the signature

    - ``determinant`` -- an integer the sign is ignored

    - ``max_scale`` -- (default: ``True``) an integer the maximum scale of a jordan block

    - ``even`` -- bool (default: ``True``)

    OUTPUT:

    A list of all global genera with the given conditions.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import all_genera_by_det
        sage: G = all_genera_by_det((4,0), 125, even=True)
        sage: G     # random
        [Genus of
        None
        Genus symbol at 2:    1^-4
        Genus symbol at 5:     1^1 5^3, Genus of
        None
        Genus symbol at 2:    1^-4
        Genus symbol at 5:     1^-2 5^1 25^-1, Genus of
        None
        Genus symbol at 2:    1^-4
        Genus symbol at 5:     1^2 5^1 25^1, Genus of
        None
        Genus symbol at 2:    1^-4
        Genus symbol at 5:     1^3 125^1]
    """
    from sage.misc.mrange import mrange_iter
    # input checks
    ZZ = IntegerRing()
    determinant = ZZ(determinant)
    sig_pair = (ZZ(sig_pair[0]), ZZ(sig_pair[1]))
    if not all([s >= 0 for s in sig_pair]):
        raise ValueError("the signature vector must be a pair of non negative integers.")
    if max_scale == None:
        max_scale = determinant
    else:
        max_scale = ZZ(max_scale)
    if type(even) != bool:
        raise ValueError("not a boolean")

    rank = sig_pair[0] + sig_pair[1]
    genera = []
    local_symbols = []
    # every global genus has a 2-adic symbol
    if determinant % 2 != 0:
        local_symbols.append(_all_p_adic_genera(2, rank, 0, 0, even=even))
    # collect the p-adic symbols
    for pn in determinant.factor():
        p = pn[0]
        det_val = pn[1]
        mscale_p = max_scale.valuation(p)
        local_symbol_p = _all_p_adic_genera(p, rank, det_val, mscale_p, even)
        local_symbols.append(local_symbol_p)
    # take the cartesian product of the collection of all possible
    # local genus symbols one for each prime
    # and check which combinations produce a global genus
    # TODO:
    # we are overcounting. Find a more
    # clever way to directly match the symbols for different primes.
    for g in mrange_iter(local_symbols):
        # create a Genus from a list of local symbols
        G = GenusSymbol_global_ring(sig_pair, g)
        # discard the empty genera
        if is_GlobalGenus(G):
            genera.append(G)
    # for testing
    genera.sort(key=lambda x: [s.symbol_tuple_list() for s in x.local_symbols()])
    return(genera)

def _all_p_adic_genera(p, rank, det_val, max_scale, even):
    r"""
    Return all `p`-adic genera with the given conditions.

    This is a helper function for :meth:`all_genera_by_det`.
    No input checks are done.

    INPUT:

    - ``p`` -- a prime number

    - ``rank`` -- the rank of this genus

    - ``det_val`` -- valuation of the determinant at p

    - ``max_scale`` -- an integer the maximal scale of a jordan block

    - ``even`` -- ``bool``; is igored if `p` is not `2`

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import _all_p_adic_genera
        sage: _all_p_adic_genera(2,3,1,2,False)
        [Genus symbol at 2:    1^-2 [2^1]_1,
         Genus symbol at 2:    1^2 [2^1]_1,
         Genus symbol at 2:    1^2 [2^1]_7,
         Genus symbol at 2:    [1^2 2^1]_3,
         Genus symbol at 2:    1^-2 [2^1]_7,
         Genus symbol at 2:    [1^-2 2^1]_7,
         Genus symbol at 2:    [1^-2 2^1]_1,
         Genus symbol at 2:    [1^2 2^1]_7,
         Genus symbol at 2:    [1^2 2^1]_5,
         Genus symbol at 2:    [1^-2 2^1]_3,
         Genus symbol at 2:    [1^-2 2^1]_5,
         Genus symbol at 2:    [1^2 2^1]_1]

    Setting a maximum scale::

        sage: _all_p_adic_genera(5, 2, 2, 1, True)
        [Genus symbol at 5:     5^-2, Genus symbol at 5:     5^2]
        sage: _all_p_adic_genera(5, 2, 2, 2, True)
        [Genus symbol at 5:     1^-1 25^-1,
         Genus symbol at 5:     1^1 25^-1,
         Genus symbol at 5:     1^-1 25^1,
         Genus symbol at 5:     1^1 25^1,
         Genus symbol at 5:     5^-2,
         Genus symbol at 5:     5^2]
    """
    from sage.misc.mrange import cantor_product
    from sage.combinat.integer_lists.invlex import IntegerListsLex
    scales_rks = [] # contains possibilities for scales and ranks
    for rkseq in IntegerListsLex(rank, length=max_scale+1):   # rank sequences
        # sum(rkseq) = rank
        # len(rkseq) = max_scale + 1
        # now assure that we get the right determinant
        d = 0
        pgensymbol = []
        for i in range(max_scale + 1):
            d += i * rkseq[i]
            # blocks of rank 0 are omitted
            if rkseq[i] != 0:
                pgensymbol.append([i, rkseq[i], 0])
        if d == det_val:
            scales_rks.append(pgensymbol)
    # add possible determinant square classes
    symbols = []
    if p != 2:
        for g in scales_rks:
            n = len(g)
            for v in cantor_product([-1, 1], repeat=n):
                g1 = deepcopy(g)
                for k in range(n):
                    g1[k][2] = v[k]
                g1 = Genus_Symbol_p_adic_ring(p, g1)
                symbols.append(g1)
    # for p == 2 we have to include determinant, even/odd, oddity
    # further restrictions apply and are defered to _blocks
    # (brute force sieving is too slow)
    # TODO: If this is too slow, enumerate only the canonical symbols.
    # as a drawback one has to reconstruct the symbol from the canonical symbol
    # this is more work for the programmer
    if p == 2:
        for g in scales_rks:
            n = len(g)
            poss_blocks = []
            for b in g:
                b += [0, 0]
                poss_blocks.append(_blocks(b, even_only=(even and b[0]==0)))
            for g1 in cantor_product(*poss_blocks):
                g1 = list(g1)
                if is_2_adic_genus(g1):
                    g1 = Genus_Symbol_p_adic_ring(p, g1)
                    # some of our symbols have the same canonical symbol
                    # thus they are equivalent - we want only one in
                    # each equivalence class
                    if not g1 in symbols:
                        symbols.append(g1)
    return(symbols)


def _blocks(b, even_only=False):
    r"""
    Return all viable `2`-adic jordan blocks with rank and scale given by ``b``

    This is a helper function for :meth:`_all_p_adic_genera`.
    It is based on the existence conditions for a modular `2`-adic genus symbol.

    INPUT:

    - ``b`` -- a list of `5` non-negative integers the first two are kept
      and all possibilities for the remaining `3` are enumerated

    - ``even_only`` -- bool (default: ``True``) if set, the blocks are even

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import _blocks
        sage: _blocks([15, 2, 0, 0, 0])
        [[15, 2, 3, 0, 0],
         [15, 2, 7, 0, 0],
         [15, 2, 1, 1, 2],
         [15, 2, 5, 1, 6],
         [15, 2, 1, 1, 6],
         [15, 2, 5, 1, 2],
         [15, 2, 7, 1, 0],
         [15, 2, 3, 1, 4]]
    """
    from copy import copy
    blocks = []
    rk = b[1]
    # recall: 2-genus_symbol is [scale, rank, det, even/odd, oddity]
    if rk == 1 and not even_only:
        for det in [1, 3, 5, 7]:
            b1 = copy(b)
            b1[2] = det
            b1[3] = 1
            b1[4] = det
            blocks.append(b1)
    elif rk == 2:
        b1 = copy(b)
        # even case
        b1[3] = 0
        b1[4] = 0
        b1[2] = 3
        blocks.append(b1)
        b1 = copy(b1)
        b1[2] = 7
        blocks.append(b1)
        # odd case
        if not even_only:
            # format (det, oddity)
            for s in [(1,2), (5,6), (1,6), (5,2), (7,0), (3,4)]:
                b1 = copy(b)
                b1[2] = s[0]
                b1[3] = 1
                b1[4] = s[1]
                blocks.append(b1)
    elif rk % 2 == 0:
        # the even case has even rank
        b1 = copy(b)
        b1[3] = 0
        b1[4] = 0
        d = (-1)**(rk//2) % 8
        for det in [d, d * (-3) % 8]:
            b1 = copy(b1)
            b1[2] = det
            blocks.append(b1)
        # odd case
        if not even_only:
            for s in [(1,2), (5,6), (1,6), (5,2), (7,0), (3,4)]:
                b1 = copy(b)
                b1[2] = s[0]*(-1)**(rk//2 -1) % 8
                b1[3] = 1
                b1[4] = s[1]
                blocks.append(b1)
            for s in [(1,4), (5,0)]:
                b1 = copy(b)
                b1[2] = s[0]*(-1)**(rk//2 - 2) % 8
                b1[3] = 1
                b1[4] = s[1]
                blocks.append(b1)
    elif rk % 2 == 1 and not even_only:
        # odd case
        for t in [1, 3, 5, 7]:
            d = (-1)**(rk//2)*t % 8
            for det in [d, -3*d % 8]:
                b1 = copy(b)
                b1[2] = det
                b1[3] = 1
                b1[4] = t
                blocks.append(b1)
    return blocks

def Genus(A, factored_determinant=None):
    r"""
    Given a nonsingular symmetric matrix `A`, return the genus of `A`.

    INPUT:

    - ``A`` -- a symmetric matrix with integer coefficients

    - ``factored_determinant`` -- (default: ``None``) a factorization object
                                  the factored determinant of ``A``

    OUTPUT:

    A :class:`GenusSymbol_global_ring` object, encoding the Conway-Sloane
    genus symbol of the quadratic form whose Gram matrix is `A`.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import Genus
        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: Genus(A)
        Genus of
        [1 1]
        [1 2]
        Genus symbol at 2:    [1^2]_2

        sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
        sage: Genus(A, A.det().factor())
        Genus of
        [2 1]
        [1 2]
        Genus symbol at 2:    1^-2
        Genus symbol at 3:     1^-1 3^-1
    """
    if factored_determinant is None:
        D = A.determinant()
        D = 2*D
        D = D.factor()
    else:
        D = factored_determinant * 2
    signature_pair = signature_pair_of_matrix(A)
    local_symbols = []
    for f in D:
        p = f[0]
        val = f[1]
        symbol = p_adic_symbol(A, p, val = val)
        G = Genus_Symbol_p_adic_ring(p, symbol)
        local_symbols.append(G)
    return GenusSymbol_global_ring(signature_pair, local_symbols, representative=A)

def LocalGenusSymbol(A, p):
    """
    Return the local symbol of `A` at the prime `p`.

    INPUT:

    - ``A`` -- a symmetric, non-singular matrix with coefficients in `\ZZ`
    - ``p`` -- an integer prime `p > 0`

    OUTPUT:

    A :class:`Genus_Symbol_p_adic_ring` object, encoding the Conway-Sloane
    genus symbol at `p` of the quadratic form whose Gram matrix is `A`.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: LocalGenusSymbol(A, 2)
        Genus symbol at 2:    [1^2]_2
        sage: LocalGenusSymbol(A, 3)
        Genus symbol at 3:     1^2

        sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
        sage: LocalGenusSymbol(A, 2)
        Genus symbol at 2:    [1^1 2^1]_2
        sage: LocalGenusSymbol(A, 3)
        Genus symbol at 3:     1^-2
    """
    val = A.determinant().valuation(p)
    symbol = p_adic_symbol(A, p, val = val)
    return Genus_Symbol_p_adic_ring(p, symbol)



def is_GlobalGenus(G):
    """
    Given a genus symbol `G` (specified by a collection of local symbols),
    return ``True`` if `G` represents the genus of a global quadratic form or lattice.

    INPUT:

    - ``G`` -- :class:`GenusSymbol_global_ring` object

    OUTPUT:

    - boolean

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import Genus, is_GlobalGenus

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G = Genus(A)
        sage: is_GlobalGenus(G)
        True

        sage: from sage.quadratic_forms.genera.genus import Genus,is_GlobalGenus
        sage: G=Genus(matrix.diagonal([2,2,2,2]))
        sage: G._local_symbols[0]._symbol=[[0,2,3,0,0],[1,2,5,1,0]]
        sage: G._representative=None
        sage: is_GlobalGenus(G)
        False

    """
    D = G.determinant()
    r, s = G.signature_pair_of_matrix()
    oddity = r - s
    for loc in G._local_symbols:
        p = loc._prime
        sym = loc._symbol
        v = sum([ss[0] * ss[1] for ss in sym])
        a = D // (p**v)
        b = Integer(prod([ss[2] for ss in sym]))
        if p == 2:
            if not is_2_adic_genus(sym):
                verbose(mesg="False in is_2_adic_genus(sym)", level=2)
                return False
            if (a*b).kronecker(p) != 1:
                verbose(mesg="False in (%s*%s).kronecker(%s)"%(a,b,p), level=2)
                return False
            oddity -= loc.excess()
        else:
            if a.kronecker(p) != b:
                verbose(mesg="False in %s.kronecker(%s) != *%s"%(a,p,b), level=2)
                return False
            oddity += loc.excess()
    if oddity % 8 != 0:
        verbose(mesg="False in oddity", level=2)
        return False
    return True



def is_2_adic_genus(genus_symbol_quintuple_list):
    """
    Given a `2`-adic local symbol (as the underlying list of quintuples)
    check whether it is the `2`-adic symbol of a `2`-adic form.

    INPUT:

    - ``genus_symbol_quintuple_list`` -- a quintuple of integers (with certain
      restrictions).

    OUTPUT:

    boolean

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2)
        sage: is_2_adic_genus(G2.symbol_tuple_list())
        True

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G3 = LocalGenusSymbol(A, 3)
        sage: is_2_adic_genus(G3.symbol_tuple_list())  ## This raises an error
        Traceback (most recent call last):
        ...
        TypeError: The genus symbols are not quintuples, so it's not a genus symbol for the prime p=2.

        sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
        sage: G2 = LocalGenusSymbol(A, 2)
        sage: is_2_adic_genus(G2.symbol_tuple_list())
        True
    """
    ## TO DO: Add explicit checking for the prime p here to ensure it's p=2... not just the quintuple checking below

    for s in genus_symbol_quintuple_list:

        ## Check that we have a quintuple (i.e. that p=2 and not p >2)
        if len(s) != 5:
            raise TypeError("The genus symbols are not quintuples, so it's not a genus symbol for the prime p=2.")

        ## Check the Conway-Sloane conditions
        if s[1] == 1:
            if s[3] == 0 or s[2] != s[4]:
                return False
        if s[1] == 2 and s[3] == 1:
            if s[2]%8 in (1,7):
               if not s[4] in (0,2,6):
                  return False
            if s[2]%8 in (3,5):
               if not s[4] in (2,4,6):
                  return False
        if (s[1] - s[4])% 2 == 1:
            return False
        if s[3] == 0 and s[4] != 0:
            return False
    return True



def canonical_2_adic_compartments(genus_symbol_quintuple_list):
    """
    Given a `2`-adic local symbol (as the underlying list of quintuples)
    this returns a list of lists of indices of the
    genus_symbol_quintuple_list which are in the same compartment.  A
    compartment is defined to be a maximal interval of Jordan
    components all (scaled) of type I (i.e. odd).

    INPUT:

    - ``genus_symbol_quintuple_list`` -- a quintuple of integers (with certain
      restrictions).

    OUTPUT:

    a list of lists of integers.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol
        sage: from sage.quadratic_forms.genera.genus import canonical_2_adic_compartments

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 1, 1, 2]]
        sage: canonical_2_adic_compartments(G2.symbol_tuple_list())
        [[0]]

        sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: canonical_2_adic_compartments(G2.symbol_tuple_list())
        [[0, 1]]

        sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
        sage: canonical_2_adic_compartments(G2.symbol_tuple_list())
        [[0, 1, 2]]

        sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 3, 0, 0]]
        sage: canonical_2_adic_compartments(G2.symbol_tuple_list())   ## No compartments here!
        []

    NOTES:

    See [Co1999]_ Conway-Sloane 3rd edition, pp. 381-382 for definitions
    and examples.
    """
    symbol = genus_symbol_quintuple_list
    compartments = []
    i = 0
    r = len(symbol)
    while i < r:
        s = symbol[i]
        if s[3] == 1:
            v = s[0]
            c = []
            while i < r and symbol[i][3] == 1 and symbol[i][0] == v:
                c.append(i)
                i += 1
                v += 1
            compartments.append(c)
        else:
            i += 1
    return compartments

def canonical_2_adic_trains(genus_symbol_quintuple_list, compartments=None):
    """
    Given a `2`-adic local symbol (as the underlying list of quintuples)
    this returns a list of lists of indices of the
    genus_symbol_quintuple_list which are in the same train.  A train
    is defined to be a maximal interval of Jordan components so that
    at least one of each adjacent pair (allowing zero-dimensional
    Jordan components) is (scaled) of type I (i.e. odd).
    Note that an interval of length one respects this condition as
    there is no pair in this interval.
    In particular, every Jordan component is part of a train.

    INPUT:

    - ``genus_symbol_quintuple_list`` -- a quintuple of integers (with certain
      restrictions).
    - ``compartments`` -- this argument is deprecated

    OUTPUT:

    a list of lists of distinct integers.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol
        sage: from sage.quadratic_forms.genera.genus import canonical_2_adic_compartments
        sage: from sage.quadratic_forms.genera.genus import canonical_2_adic_trains

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 1, 1, 2]]
        sage: canonical_2_adic_trains(G2.symbol_tuple_list())
        [[0]]

        sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: canonical_2_adic_compartments(G2.symbol_tuple_list())
        [[0, 1]]

        sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
        sage: canonical_2_adic_trains(G2.symbol_tuple_list())
        [[0, 1, 2]]

        sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 3, 0, 0]]
        sage: canonical_2_adic_trains(G2.symbol_tuple_list())
        [[0]]
        sage: symbol = [[0, 1,  1, 1, 1],[1, 2, -1, 0, 0],[2, 1,  1, 1, 1],[3, 1,  1, 1, 1],[4, 1,  1, 1, 1],[5, 2, -1, 0, 0],[7, 1,  1, 1, 1],[10, 1, 1, 1, 1],[11, 1, 1, 1, 1],[12, 1, 1, 1, 1]]
        sage: canonical_2_adic_trains(symbol)
        [[0, 1, 2, 3, 4, 5], [6], [7, 8, 9]]

    Check that :trac:`24818` is fixed::

        sage: symbol = [[0, 1,  1, 1, 1],[1, 3, 1, 1, 1]]
        sage: canonical_2_adic_trains(symbol)
        [[0, 1]]

    .. NOTE::

        See [Co1999]_, pp. 381-382 for definitions and examples.

    """
    if compartments is not None:
        from sage.misc.superseded import deprecation
        deprecation(23955, "the compartments keyword has been deprecated")

    # avoid a special case for the end of symbol
    # if a jordan component has rank zero it is considered even.
    symbol = genus_symbol_quintuple_list
    symbol.append([symbol[-1][0]+1, 0, 1, 0, 0]) #We have just modified the input globally!
    # Hence, we have to remove the last entry of symbol at the end.
    try:

        trains = []
        new_train = [0]
        for i in range(1,len(symbol)-1):
            # start a new train if there are two adjacent even symbols
            prev, cur = symbol[i-1:i+1]
            if  cur[0] - prev[0] > 2:
                trains.append(new_train)
                new_train = [i]    # create a new train starting at
            elif (cur[0] - prev[0] == 2) and cur[3]*prev[3] == 0:
                trains.append(new_train)
                new_train = [i]
            elif prev[3] == 0 and cur[3] == 0:
                trains.append(new_train)
                new_train = [i]
            else:
                # there is an odd jordan block adjacent to this jordan block
                # the train continues
                new_train.append(i)
        # the last train was never added.
        trains.append(new_train)
        return trains
    finally:
        #revert the input list to its original state
        symbol.pop()

def canonical_2_adic_reduction(genus_symbol_quintuple_list):
    """
    Given a `2`-adic local symbol (as the underlying list of quintuples)
    this returns a canonical `2`-adic symbol (again as a raw list of
    quintuples of integers) which has at most one minus sign per train
    and this sign appears on the smallest dimensional Jordan component
    in each train.  This results from applying the "sign-walking" and
    "oddity fusion" equivalences.

    INPUT:

    - ``genus_symbol_quintuple_list`` -- a quintuple of integers (with certain
      restrictions)

    - ``compartments`` -- a list of lists of distinct integers (optional)

    OUTPUT:

    a list of lists of distinct integers.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol
        sage: from sage.quadratic_forms.genera.genus import canonical_2_adic_reduction

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 1, 1, 2]]
        sage: canonical_2_adic_reduction(G2.symbol_tuple_list())
        [[0, 2, 1, 1, 2]]

        sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: canonical_2_adic_reduction(G2.symbol_tuple_list())   ## Oddity fusion occurred here!
        [[0, 1, 1, 1, 2], [1, 1, 1, 1, 0]]

        sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
        sage: canonical_2_adic_reduction(G2.symbol_tuple_list())   ## Oddity fusion occurred here!
        [[1, 2, -1, 1, 6], [2, 1, 1, 1, 0], [3, 1, 1, 1, 0]]

        sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
        sage: G2 = LocalGenusSymbol(A, 2); G2.symbol_tuple_list()
        [[0, 2, 3, 0, 0]]
        sage: canonical_2_adic_reduction(G2.symbol_tuple_list())
        [[0, 2, -1, 0, 0]]

    .. NOTE::

        See [Co1999]_ Conway-Sloane 3rd edition, pp. 381-382 for definitions and examples.

    .. TODO::

        Add an example where sign walking occurs!
    """
    # Protect the input from unwanted modification
    genus_symbol_quintuple_list = deepcopy(genus_symbol_quintuple_list)
    canonical_symbol = genus_symbol_quintuple_list
    # Canonical determinants:
    for i in range(len(genus_symbol_quintuple_list)):
        d = genus_symbol_quintuple_list[i][2]
        if d in (1,7):
            canonical_symbol[i][2] = 1
        else:
            canonical_symbol[i][2] = -1
    # Oddity fusion:
    compartments = canonical_2_adic_compartments(genus_symbol_quintuple_list)
    for compart in compartments:
        oddity = sum([ genus_symbol_quintuple_list[i][4] for i in compart ]) % 8
        for i in compart:
            genus_symbol_quintuple_list[i][4] = 0
        genus_symbol_quintuple_list[compart[0]][4] = oddity
    verbose(mesg="End oddity fusion: %s" %canonical_symbol, level=2)
    # Sign walking:
    trains = canonical_2_adic_trains(genus_symbol_quintuple_list)
    for train in trains:
        t = len(train)
        for i in range(t-1):
            t1 = train[t-i-1]
            if canonical_symbol[t1][2] == -1:
                canonical_symbol[t1][2] = 1
                canonical_symbol[t1-1][2] *= -1
                for compart in compartments:
                    if t1-1 in compart or t1 in compart:
                        o = canonical_symbol[compart[0]][4]
                        canonical_symbol[compart[0]][4] = (o+4) % 8
    verbose(mesg="End sign walking: %s" %canonical_symbol, level=2)
    return canonical_symbol


def basis_complement(B):
    """
    Given an echelonized basis matrix `B` (over a field), calculate a
    matrix whose rows form a basis complement (to the rows of `B`).

    INPUT:

    - ``B`` -- matrix over a field in row echelon form

    OUTPUT:

    a rectangular matrix over a field

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import basis_complement

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,1])
        sage: B = A.kernel().echelonized_basis_matrix(); B
        [ 1 -1]
        sage: basis_complement(B)
        [0 1]
    """
    F = B.parent().base_ring()
    m = B.nrows()
    n = B.ncols()
    C = MatrixSpace(F,n-m,n,sparse=True)(0)
    k = 0
    l = 0
    for i in range(m):
        for j in range(k,n):
             if B[i,j] == 0:
                 C[l,j] = 1
                 l += 1
             else:
                 k = j+1
                 break
    for j in range(k,n):
        C[l+j-k,j] = 1
    return C



def signature_pair_of_matrix(A):
    """
    Computes the signature pair `(p, n)` of a non-degenerate symmetric
    matrix, where

    - `p` is the number of positive eigenvalues of `A`
    - `n` is the number of negative eigenvalues of `A`

    INPUT:

    - ``A`` -- symmetric matrix (assumed to be non-degenerate)

    OUTPUT:

    - `(p, n)` -- a pair (tuple) of integers.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import signature_pair_of_matrix

        sage: A = Matrix(ZZ, 2, 2, [-1,0,0,3])
        sage: signature_pair_of_matrix(A)
        (1, 1)

        sage: A = Matrix(ZZ, 2, 2, [-1,1,1,7])
        sage: signature_pair_of_matrix(A)
        (1, 1)

        sage: A = Matrix(ZZ, 2, 2, [3,1,1,7])
        sage: signature_pair_of_matrix(A)
        (2, 0)

        sage: A = Matrix(ZZ, 2, 2, [-3,1,1,-11])
        sage: signature_pair_of_matrix(A)
        (0, 2)


        sage: A = Matrix(ZZ, 2, 2, [1,1,1,1])
        sage: signature_pair_of_matrix(A)
        Traceback (most recent call last):
        ...
        ArithmeticError: given matrix is not invertible
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm
    s_vec = QuadraticForm(A.base_extend(A.base_ring().fraction_field())).signature_vector()

    # Check that the matrix is non-degenerate (i.e. no zero eigenvalues)
    if s_vec[2]:
        raise ArithmeticError("given matrix is not invertible")

    # Return the pair (p,n)
    return s_vec[:2]


def p_adic_symbol(A, p, val):
    """
    Given a symmetric matrix `A` and prime `p`, return the genus symbol at `p`.

    .. TODO::

        Some description of the definition of the genus symbol.

    INPUT:

    - ``A`` -- symmetric matrix with integer coefficients
    - ``p`` -- prime number > 0
    - ``val`` -- integer >= 0; valuation of the maximal elementary divisor of `A`
      needed to obtain enough precision. Calculation is modulo `p` to the ``val+3``.

    OUTPUT:

    a list of lists of integers

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import p_adic_symbol

        sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
        sage: p_adic_symbol(A, 2, 2)
        [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]

        sage: p_adic_symbol(A, 3, 1)
        [[0, 3, 1], [1, 1, -1]]

    """
    if p % 2 == 0:
        return two_adic_symbol(A, val)
    m0 = min([ c.valuation(p) for c in A.list() ])
    q = p**m0
    n = A.nrows()
    A = MatrixSpace(IntegerRing(),n,n)([ c // q for c in A.list() ])
    A_p = MatrixSpace(FiniteField(p),n,n)(A)
    B_p = A_p.kernel().echelonized_basis_matrix()
    if B_p.nrows() == 0:
        e0 = Integer(A_p.det()).kronecker(p)
        n0 = A.nrows()
        return [ [m0,n0,e0] ]
    else:
        C_p = basis_complement(B_p)
        e0 = Integer((C_p*A_p*C_p.transpose()).det()).kronecker(p)
        n0 = C_p.nrows()
        sym = [ [0,n0,e0] ]
    r = B_p.nrows()
    B = MatrixSpace(IntegerRing(),r,n)(B_p)
    C = MatrixSpace(IntegerRing(),n-r,n)(C_p)
    # Construct the blocks for the Jordan decomposition [F,X;X,A_new]
    F = MatrixSpace(RationalField(),n-r,n-r)(C*A*C.transpose())
    U = F**-1
    d = LCM([ c.denominator() for c in U.list() ])
    R = IntegerRing().quotient_ring(Integer(p)**(val+3))
    u = R(d)**-1
    MatR = MatrixSpace(R,n-r,n-r)
    MatZ = MatrixSpace(IntegerRing(),n-r,n-r)
    U = MatZ(MatR(MatZ(U*d))*u)
    # X = C*A*B.transpose()
    # A = B*A*B.transpose() - X.transpose()*U*X
    X = C*A
    A = B*(A - X.transpose()*U*X)*B.transpose()
    return [ [s[0]+m0] + s[1:] for s in sym + p_adic_symbol(A, p, val) ]



def is_even_matrix(A):
    """
    Determines if the integral symmetric matrix `A` is even
    (i.e. represents only even numbers).  If not, then it returns the
    index of an odd diagonal entry.  If it is even, then we return the
    index -1.

    INPUT:

    - ``A`` -- symmetric integer matrix

    OUTPUT:

    a pair of the form (boolean, integer)

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import is_even_matrix

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,1])
        sage: is_even_matrix(A)
        (False, 0)

        sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
        sage: is_even_matrix(A)
        (True, -1)
    """
    for i in range(A.nrows()):
        if A[i,i]%2 == 1:
            return False, i
    return True, -1



def split_odd(A):
    """
    Given a non-degenerate Gram matrix `A (\mod 8)`, return a splitting
    ``[u] + B`` such that u is odd and `B` is not even.

    INPUT:

    - ``A`` -- an odd symmetric matrix with integer coefficients (which admits a
      splitting as above).

    OUTPUT:

    a pair ``(u, B)`` consisting of an odd integer `u` and an odd
    integral symmetric matrix `B`.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import is_even_matrix
        sage: from sage.quadratic_forms.genera.genus import split_odd

        sage: A = Matrix(ZZ, 2, 2, [1,2,2,3])
        sage: is_even_matrix(A)
        (False, 0)
        sage: split_odd(A)
        (1, [-1])

        sage: A = Matrix(ZZ, 2, 2, [1,2,2,5])
        sage: split_odd(A)
        (1, [1])

        sage: A = Matrix(ZZ, 2, 2, [1,1,1,1])
        sage: is_even_matrix(A)
        (False, 0)
        sage: split_odd(A)      ## This fails because no such splitting exists. =(
        Traceback (most recent call last):
        ...
        RuntimeError: The matrix A does not admit a non-even splitting.

        sage: A = Matrix(ZZ, 2, 2, [1,2,2,6])
        sage: split_odd(A)      ## This fails because no such splitting exists. =(
        Traceback (most recent call last):
        ...
        RuntimeError: The matrix A does not admit a non-even splitting.

    """
    n0 = A.nrows()
    if n0 == 1:
       return A[0,0], MatrixSpace(IntegerRing(),0,A.ncols())([])
    even, i = is_even_matrix(A)
    R = A.parent().base_ring()
    C = MatrixSpace(R,n0-1,n0)(0)
    u = A[i,i]
    for j in range(n0-1):
        if j < i:
            C[j,j] = 1
            C[j,i] = -A[j,i]*u
        else:
            C[j,j+1] = 1
            C[j,i] = -A[j+1,i]*u
        B = C*A*C.transpose()
    even, j = is_even_matrix(B)
    if even:
        I = A.parent()(1)
        # TODO: we could manually (re)construct the kernel here...
        if i == 0:
            I[1,0] = 1 - A[1,0]*u
            i = 1
        else:
            I[0,i] = 1 - A[0,i]*u
            i = 0
        A = I*A*I.transpose()
        u = A[i,i]
        C = MatrixSpace(R,n0-1,n0)(0)
        for j in range(n0-1):
            if j < i:
               C[j,j] = 1
               C[j,i] = -A[j,i]*u
            else:
                C[j,j+1] = 1
                C[j,i] = -A[j+1,i]*u
            B = C*A*C.transpose()
    even, j = is_even_matrix(B)
    if even:
        print("B:")
        print(B)
        raise RuntimeError("The matrix A does not admit a non-even splitting.")
    return u, B



def trace_diag_mod_8(A):
    """
    Return the trace of the diagonalised form of `A` of an integral
    symmetric matrix which is diagonalizable `\mod 8`.  (Note that since
    the Jordan decomposition into blocks of size `<=` 2 is not unique
    here, this is not the same as saying that `A` is always diagonal in
    any `2`-adic Jordan decomposition!)

    INPUT:

    - ``A`` -- symmetric matrix with coefficients in `\ZZ` which is odd in
      `\ZZ/2\ZZ` and has determinant not divisible by `8`.

    OUTPUT:

    an integer

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import is_even_matrix
        sage: from sage.quadratic_forms.genera.genus import split_odd
        sage: from sage.quadratic_forms.genera.genus import trace_diag_mod_8

        sage: A = Matrix(ZZ, 2, 2, [1,2,2,3])
        sage: is_even_matrix(A)
        (False, 0)
        sage: split_odd(A)
        (1, [-1])
        sage: trace_diag_mod_8(A)
        0

        sage: A = Matrix(ZZ, 2, 2, [1,2,2,5])
        sage: split_odd(A)
        (1, [1])
        sage: trace_diag_mod_8(A)
        2
    """
    tr = 0
    while A.nrows() > 0:
       u, A = split_odd(A)
       tr += u
    return IntegerRing()(tr)



def two_adic_symbol(A, val):
    """
    Given a symmetric matrix `A` and prime `p`, return the genus symbol at `p`.

    The genus symbol of a component 2^m*f is of the form ``(m,n,s,d[,o])``,
    where

    - m = valuation of the component
    - n = dimension of f
    - d = det(f) in {1,3,5,7}
    - s = 0 (or 1) if even (or odd)
    - o = oddity of f (= 0 if s = 0) in `Z/8Z`

    INPUT:

    - ``A`` -- symmetric matrix with integer coefficients, non-degenerate
    - ``val`` -- integer >=0; valuation of maximal 2-elementary divisor

    OUTPUT:

    a list of lists of integers (representing a Conway-Sloane `2`-adic symbol)

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import two_adic_symbol

        sage: A = diagonal_matrix(ZZ, [1,2,3,4])
        sage: two_adic_symbol(A, 2)
        [[0, 2, 3, 1, 4], [1, 1, 1, 1, 1], [2, 1, 1, 1, 1]]

    """
    m0 = min([ c.valuation(2) for c in A.list() ])
    q = 2**m0
    A = A.parent()([ c // q for c in A.list() ])
    ZZ = IntegerRing()
    n = A.nrows()
    A_2 = MatrixSpace(FiniteField(2),n,n)(A)
    K_2 = A_2.kernel()
    R_8 = ZZ.quotient_ring(Integer(8))

    ## Deal with the matrix being non-degenerate mod 2.
    if K_2.dimension() == 0:
        A_8 = MatrixSpace(R_8,n)(A)
        n0 = A.nrows()
        # d0 = ZZ(A_8.determinant()) # no determinant over Z/8Z
        d0 = ZZ(R_8(MatrixSpace(ZZ,n)(A_8).determinant()))
        if d0 == 0:    ## SANITY CHECK: The mod 8 determinant shouldn't be zero.
            print("A:")
            print(A)
            assert False
        even, i = is_even_matrix(A_2)    ## Determine whether the matrix is even or odd.
        if even:
            return [ [m0,n0,d0,0,0] ]
        else:
            tr8 = trace_diag_mod_8(A_8)  ## Here we already know that A_8 is odd and diagonalizable mod 8.
            return [ [m0,n0,d0,1,tr8] ]

    ## Deal with the matrix being degenerate mod 2.
    else:
        B_2 = K_2.echelonized_basis_matrix()
        C_2 = basis_complement(B_2)
        n0 = C_2.nrows()
        C = MatrixSpace(ZZ,n0,n)(C_2)
        A_new = C*A*C.transpose()
        # compute oddity modulo 8:
        A_8 = MatrixSpace(R_8,n0,n0)(A_new)
        # d0 = A_8.det() # no determinant over Z/8Z
        d0 = ZZ(R_8(MatrixSpace(ZZ,n0,n0)(A_8).determinant()))
        if d0 == 0:
            print("A:")
            print(A_new)
            assert False
        even, i = is_even_matrix(A_new)
        if even:
            sym = [ [0,n0,d0,0,0] ]
        else:
            tr8 = trace_diag_mod_8(A_8)
            sym = [ [0,n0,d0,1,tr8] ]
    r = B_2.nrows()
    B = MatrixSpace(ZZ,r,n)(B_2)
    C = MatrixSpace(IntegerRing(),n-r,n)(C_2)
    F = MatrixSpace(RationalField(),n-r,n-r)(C*A*C.transpose())
    U = F**-1
    d = LCM([ c.denominator() for c in U.list() ])
    R = IntegerRing().quotient_ring(Integer(2)**(val+3))
    u = R(d)**-1
    MatR = MatrixSpace(R,n-r,n-r)
    MatZ = MatrixSpace(IntegerRing(),n-r,n-r)
    U = MatZ(MatR(MatZ(U*d))*u)
    X = C*A
    A = B*(A - X.transpose()*U*X)*B.transpose()
    return [ [s[0]+m0] + s[1:] for s in sym + two_adic_symbol(A, val) ]


class Genus_Symbol_p_adic_ring(object):
    r"""
    Local genus symbol over a p-adic ring.

    The genus symbol of a component `p^m A` for odd prime `= p` is of the
    form `(m,n,d)`, where

    - `m` = valuation of the component
    - `n` = rank of A
    - `d = det(A) \in \{1,u\}` for a normalized quadratic non-residue `u`.

    The genus symbol of a component `2^m A` is of the form `(m, n, s, d, o)`,
    where

    - `m` = valuation of the component
    - `n` = rank of `A`
    - `d` = det(A) in `\{1,3,5,7\}`
    - `s` = 0 (or 1) if even (or odd)
    - `o` = oddity of `A` (= 0 if s = 0) in `Z/8Z`
          = the trace of the diagonalization of `A`

    The genus symbol is a list of such symbols (ordered by `m`) for each
    of the Jordan blocks `A_1,...,A_t`.

    Reference: [Co1999]_ Conway and Sloane 3rd edition, Chapter 15, Section 7.


    WARNING/NOTE: This normalization seems non-standard, and we
    should review this entire class to make sure that we have our
    doubling conventions straight throughout!  This is especially
    noticeable in the determinant and excess methods!!

    INPUT:

    - ``prime`` -- a prime integer `> 0`
    - ``symbol`` -- the list of invariants for Jordan blocks `A_t,...,A_t` given
      as a list of lists of integers

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
        sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

        sage: A = diagonal_matrix(ZZ, [1,2,3,4])
        sage: p = 2
        sage: s2 = p_adic_symbol(A, p, 2); s2
        [[0, 2, 3, 1, 4], [1, 1, 1, 1, 1], [2, 1, 1, 1, 1]]
        sage: G2 = Genus_Symbol_p_adic_ring(p,s2);G2
        Genus symbol at 2:    [1^-2 2^1 4^1]_6

        sage: A = diagonal_matrix(ZZ, [1,2,3,4])
        sage: p = 3
        sage: s3 = p_adic_symbol(A, p, 1); s3
        [[0, 3, -1], [1, 1, 1]]
        sage: G3 = Genus_Symbol_p_adic_ring(p,s3);G3
        Genus symbol at 3:     1^-3 3^1
    """
    def __init__(self, prime, symbol, check = True):
        """
        Create the local genus symbol of given prime and local invariants.


        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
            sage: A = diagonal_matrix(ZZ, [1,2,3,4])
            sage: p = 2
            sage: s2 = p_adic_symbol(A, p, 2); s2
            [[0, 2, 3, 1, 4], [1, 1, 1, 1, 1], [2, 1, 1, 1, 1]]
            sage: G2 = Genus_Symbol_p_adic_ring(p,s2);G2
            Genus symbol at 2:    [1^-2 2^1 4^1]_6

            sage: A = diagonal_matrix(ZZ, [1,2,3,4])
            sage: p = 3
            sage: s3 = p_adic_symbol(A, p, 1); s3
            [[0, 3, -1], [1, 1, 1]]
            sage: G3 = Genus_Symbol_p_adic_ring(p,s3);G3
            Genus symbol at 3:     1^-3 3^1
            sage: G2 == loads(dumps(G2))
            True
            sage: G3 == loads(dumps(G3))
            True
        """
        if check:
           pass
        self._prime = ZZ(prime)
        self._symbol = symbol
        self._canonical_symbol = None

    def __repr__(self):
        r"""
        String representation for the `p`-adic genus symbol

        OUTPUT:

        a string

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
            sage: symbol = [[0, 4, -1, 0, 0],[1, 2, 1, 1, 2],[2, 1, 1, 1, 1],[4, 4, 1, 0, 0],[5, 1, 1, 1, 1]]
            sage: g = Genus_Symbol_p_adic_ring(2,symbol)
            sage: g
            Genus symbol at 2:    1^-4 [2^2 4^1]_3:16^4 [32^1]_1

        TESTS:

        Check that :trac:`25776` is fixed::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: G = Genus(matrix.diagonal([2,2,64]))
            sage: G
            Genus of
            [ 2  0  0]
            [ 0  2  0]
            [ 0  0 64]
            Genus symbol at 2:    [2^2]_2:[64^1]_1
        """
        p=self._prime
        CS_string = ""
        if p==2:
            CS = self.canonical_symbol()
            for train in self.trains():
                #mark the beginning of a train with a colon
                CS_string += " :"
                #collect the indices where compartments begin and end
                compartment_begins = []
                compartment_ends = []
                for comp in self.compartments():
                    compartment_begins.append(comp[0])
                    compartment_ends.append(comp[-1])

                for block_index in train:
                    if block_index in compartment_begins:
                        #mark the beginning of this compartment with [
                        CS_string += "["
                    block = CS[block_index]
                    block_string = "%s^%s " % (p**block[0],block[2]*block[1])
                    CS_string += block_string
                    if block_index in compartment_ends:
                        #close this compartment with ] and remove a space
                        CS_string = CS_string[:-1] + "]"
                        # the oddity belongs to the compartment
                        # and is saved in its first block
                        i = compartment_ends.index(block_index)
                        compartment_start = compartment_begins[i]
                        oddity = CS[compartment_start][4]
                        CS_string +="_%s" % oddity
            # remove the first colon
            CS_string = CS_string[2:]
            # remove some unnecessary whitespace
            CS_string = CS_string.replace(" :",":")

        else:
            for s in self._symbol:
                CS_string += " %s^%s" % (p**s[0], s[2]*s[1])
        return "Genus symbol at %s:    %s" % (p, CS_string)

    def _latex_(self):
        """
        The LaTeX representation of this local genus symbol.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
            sage: symbol = [[0, 4, -1, 0, 0],[1, 2, 1, 1, 2],[2, 1, 1, 1, 1],[4, 4, 1, 0, 0],[5, 1, 1, 1, 1]]
            sage: g = Genus_Symbol_p_adic_ring(2,symbol)
            sage: g._canonical_symbol = [[0, 4, 1, 0, 0],[1, 2, 1, 1, 3],[2, 1, 1, 1, 0],[4, 4, 1, 0, 0],[5, 1, 1, 1, 1]]
            sage: g._latex_()
            '\\mbox{Genus symbol at } 2\\mbox{: }1^{4} [2^{2} 4^{1}]_{3} :16^{4} [32^{1}]_{1}'
        """
        p=self._prime
        CS_string = ""
        if p==2:
            CS = self.canonical_symbol()
            for train in self.trains():
                # mark the beginning of a train with a colon
                CS_string += " :"
                # collect the indices where compartments begin and end
                compartment_begins = []
                compartment_ends = []
                for comp in self.compartments():
                    compartment_begins.append(comp[0])
                    compartment_ends.append(comp[-1])

                for block_index in train:
                    if block_index in compartment_begins:
                        # mark the beginning of this compartment with [
                        CS_string += "["
                    block = CS[block_index]
                    block_string = "%s^{%s} " % (p**block[0],block[2]*block[1])
                    CS_string += block_string
                    if block_index in compartment_ends:
                        # close this compartment with ] and remove a space
                        CS_string = CS_string[:-1] + "]"
                        # the oddity belongs to the compartment
                        # and is saved in its first block
                        i = compartment_ends.index(block_index)
                        compartment_start = compartment_begins[i]
                        oddity = CS[compartment_start][4]
                        CS_string +="_{%s}" % oddity
            #remove the first colon
            CS_string = CS_string[2:]

        else:
            for s in self._symbol:
                CS_string += " {%s}^{%s}" % (p**s[0], s[2]*s[1])
        return "\\mbox{Genus symbol at } %s\\mbox{: }%s" % (p,CS_string)

    def __eq__(self, other):
        """
        Determines if two genus symbols are equal (not just equivalent!).

        INPUT:

        - other -- a :class:`Genus_Symbol_p_adic_ring` object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = diagonal_matrix(ZZ, [1,2,3,4])
            sage: p = 2
            sage: G2 =  Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2))
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 1))

            sage: G2 == G3
            False
            sage: G3 == G2
            False
            sage: G2 == G2
            True
            sage: G3 == G3
            True

        """
        p = self._prime
        if p != other._prime:
            return False
        return self.canonical_symbol() == other.canonical_symbol()


    def __ne__(self, other):
        """
        Determines if two genus symbols are unequal (not just inequivalent!).

        INPUT:

        - other -- a :class:`Genus_Symbol_p_adic_ring` object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = diagonal_matrix(ZZ, [1,2,3,4])
            sage: p = 2
            sage: G2 =  Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2))
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 1))

            sage: G2 != G3
            True
            sage: G3 != G2
            True
            sage: G2 != G2
            False
            sage: G3 != G3
            False

        """
        return not self == other


    ## Added these two methods to make this class iterable...
    #def  __getitem__(self, i):
    #    return self._symbol[i]
    #
    #def len(self):
    #    return len(self._symbol)
    ## ------------------------------------------------------

    def automorphous_numbers(self):
        r"""
        Return generators of the automorphous square classes at this prime.

        A `p`-adic square class `r` is called automorphous if it is
        the spinor norm of a proper `p`-adic integral automorphism of this form.
        These classes form a group. See [CS]_ chapter 15 §9.6 for details.

        OUTPUT:

        - a list of integers representing the square classes of generators of
          the automorphous numbers

        EXAMPLES:

        The following examples are given in
        [Co1999]_ 3rd edition, Chapter 15, 9.6 pp. 392::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: A = matrix.diagonal([3,16])
            sage: G = Genus(A)
            sage: sym2 = G.local_symbols()[0]
            sage: sym2
            Genus symbol at 2:    [1^-1]_3:[16^1]_1
            sage: sym2.automorphous_numbers()
            [3, 5]

            sage: A = matrix(ZZ,3,[2,1,0, 1,2,0, 0,0,18])
            sage: G = Genus(A)
            sage: sym = G.local_symbols()
            sage: sym[0]
            Genus symbol at 2:    1^-2 [2^1]_1
            sage: sym[0].automorphous_numbers()
            [1, 3, 5, 7]
            sage: sym[1]
            Genus symbol at 3:     1^-1 3^-1 9^-1
            sage: sym[1].automorphous_numbers()
            [1, 3]

        Note that the generating set given is not minimal.
        The first supplementation rule is used here::

            sage: A = matrix.diagonal([2,2,4])
            sage: G = Genus(A)
            sage: sym = G.local_symbols()
            sage: sym[0]
            Genus symbol at 2:    [2^2 4^1]_3
            sage: sym[0].automorphous_numbers()
            [1, 2, 3, 5, 7]

        but not there::

            sage: A = matrix.diagonal([2,2,32])
            sage: G = Genus(A)
            sage: sym = G.local_symbols()
            sage: sym[0]
            Genus symbol at 2:    [2^2]_2:[32^1]_1
            sage: sym[0].automorphous_numbers()
            [1, 2, 5]

        Here the second supplementation rule is used::

            sage: A = matrix.diagonal([2,2,64])
            sage: G = Genus(A)
            sage: sym = G.local_symbols()
            sage: sym[0]
            Genus symbol at 2:    [2^2]_2:[64^1]_1
            sage: sym[0].automorphous_numbers()
            [1, 2, 5]
        """
        from .normal_form import collect_small_blocks, _min_nonsquare
        automorphs = []
        sym = self.symbol_tuple_list()
        G = self.gram_matrix().change_ring(ZZ)
        p = self.prime()
        if p != 2:
            up = ZZ(_min_nonsquare(p))
            I = G.diagonal()
            for r in I:
                # We need to consider all pairs in I
                # since at most 2 elements are part of a pair
                # we need need at most 2 of each type
                if I.count(r) > 2:
                    I.remove(r)
            # products of all pairs
            for r1 in I:
                for r2 in I:
                    automorphs.append(r1*r2)
            # supplement (i)
            for block in sym:
                if block[1] >= 2:
                    automorphs.append(up)
                    break
            # normalize the square classes and remove duplicates
            automorphs1 = set()
            for s in automorphs:
                u = 1
                if s.prime_to_m_part(p).kronecker(p) == -1:
                    u = up
                v = (s.valuation(p) % 2)
                sq = u * p**v
                automorphs1.add(sq)
            return list(automorphs1)

        # p = 2
        I = []
        II = []
        for block in collect_small_blocks(G):
            if block.ncols() == 1:
                u = block[0,0]
                if I.count(u) < 2:
                    I.append(block[0,0])
            else: # rank2
                q = block[0,1]
                II += [2*q, 3*2*q, 5*2*q, 7*2*q]

        L = I + II
        # We need to consider all pairs in L
        # since at most 2 elements are part of a pair
        # we need need at most 2 of each type
        for r in L:     # remove triplicates
            if L.count(r) > 2:
                L.remove(r)
        n = len(L)
        for i in range(n):
            for j in range(i):
                r = L[i]*L[j]
                automorphs.append(r)

        # supplement (i)
        for k in range(len(sym)):
            s = sym[k:k+3]
            if sum([b[1] for b in s if b[0] - s[0][0] < 4]) >= 3:
                automorphs += [ZZ(1), ZZ(3), ZZ(5), ZZ(7)]
            break

        # supplement (ii)
        I.sort(key=lambda x: x.valuation(2))
        n = len(I)
        for i in range(n):
            for j in range(i):
                r = I[i] / I[j]
                v, u = r.val_unit(ZZ(2))
                u = u % 8
                assert v >= 0
                if v==0 and u==1:
                    automorphs.append(ZZ(2))
                if v==0 and u==5:
                    automorphs.append(ZZ(6))
                if v in [0, 2, 4]:  # this overlaps with the first two cases!
                    automorphs.append(ZZ(5))
                if v in [1, 3] and u in [1, 5]:
                    automorphs.append(ZZ(3))
                if v in [1, 3] and u in [3, 7]:
                    automorphs.append(ZZ(7))

        # normalize the square classes and remove duplicates
        automorphs1 = set()
        for s in automorphs:
            v, u = s.val_unit(ZZ(2))
            v = v % 2
            u = u % 8
            sq = u * 2**v
            automorphs1.add(sq)
        return list(automorphs1)

    def canonical_symbol(self):
        """
        Return (and cache) the canonical p-adic genus symbol.  This is
        only really affects the `2`-adic symbol, since when `p > 2` the
        symbol is already canonical.

        OUTPUT:

        a list of lists of integers

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = Matrix(ZZ, 2, 2, [1,1,1,2])
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2.symbol_tuple_list()
            [[0, 2, 1, 1, 2]]
            sage: G2.canonical_symbol()
            [[0, 2, 1, 1, 2]]

            sage: A = Matrix(ZZ, 2, 2, [1,0,0,2])
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2.symbol_tuple_list()
            [[0, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: G2.canonical_symbol()   ## Oddity fusion occurred here!
            [[0, 1, 1, 1, 2], [1, 1, 1, 1, 0]]

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2.symbol_tuple_list()
            [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
            sage: G2.canonical_symbol()   ## Oddity fusion occurred here!
            [[1, 2, -1, 1, 6], [2, 1, 1, 1, 0], [3, 1, 1, 1, 0]]

            sage: A = Matrix(ZZ, 2, 2, [2,1,1,2])
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2.symbol_tuple_list()
            [[0, 2, 3, 0, 0]]
            sage: G2.canonical_symbol()
            [[0, 2, -1, 0, 0]]


            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G3.symbol_tuple_list()
            [[0, 3, 1], [1, 1, -1]]
            sage: G3.canonical_symbol()
            [[0, 3, 1], [1, 1, -1]]

        .. NOTE::

            See [Co1999]_ Conway-Sloane 3rd edition, pp. 381-382 for definitions and examples.

        .. TODO::

            Add an example where sign walking occurs!
        """
        symbol = self._symbol
        if self._prime == 2:
            if self._canonical_symbol is None:
                self._canonical_symbol = canonical_2_adic_reduction(symbol)
            return self._canonical_symbol
        else:
            return self._symbol


    def gram_matrix(self, check=True):
        r"""
        Return a gram matrix of a representative of this local genus.

        INPUT:

        - check (default: ``True``) -- double check the result

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2))
            sage: G2.gram_matrix()
            [2 0|0|0]
            [0 6|0|0]
            [---+-+-]
            [0 0|4|0]
            [---+-+-]
            [0 0|0|8]
        """
        G = []
        p = self._prime
        symbol = self.symbol_tuple_list()
        for block in symbol:
            G.append(_gram_from_jordan_block(p, block))
        G = matrix.block_diagonal(G)
        # check calculation
        if check:
            symG = p_adic_symbol(G, p, symbol[-1][0])
            assert Genus_Symbol_p_adic_ring(p, symG) == self, "oops"
        return G.change_ring(ZZ)

    def prime(self):
        r"""
        Return the prime number `p` of this `p`-adic local symbol.

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol
            sage: M1 = matrix(ZZ,[2])
            sage: p = 2
            sage: G0 = LocalGenusSymbol(M1, 2)
            sage: G0.prime()
            2
        """
        return self._prime

    def is_even(self):
        r"""
        Return if the underlying `p`-adic lattice is even.

        If `p` is odd, every lattice is even.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import LocalGenusSymbol
            sage: M0 = matrix(ZZ,[1])
            sage: G0 = LocalGenusSymbol(M0, 2)
            sage: G0.is_even()
            False
            sage: G1 = LocalGenusSymbol(M0, 3)
            sage: G1.is_even()
            True
            sage: M2 = matrix(ZZ,[2])
            sage: G2 = LocalGenusSymbol(M2, 2)
            sage: G2.is_even()
            True
        """
        if self.prime() != 2:
            return True
        sym = self.symbol_tuple_list()[0]
        return sym[0] > 0 or sym[3]==0

    def symbol_tuple_list(self):
        """
        Return a copy of the underlying list of lists of integers
        defining the genus symbol.

        OUTPUT:

        a list of lists of integers

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G3
            Genus symbol at 3:     1^3 3^-1
            sage: G3.symbol_tuple_list()
            [[0, 3, 1], [1, 1, -1]]
            sage: type(G3.symbol_tuple_list())
            <... 'list'>

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            sage: G2.symbol_tuple_list()
            [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
            sage: type(G2.symbol_tuple_list())
            <... 'list'>

        """
        return deepcopy(self._symbol)

    def number_of_blocks(self):
        """
        Return the number of positive dimensional symbols/Jordan blocks.

        OUTPUT:

        integer >= 0

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2.symbol_tuple_list()
            [[1, 2, 3, 1, 4], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1]]
            sage: G2.number_of_blocks()
            3

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G3.symbol_tuple_list()
            [[0, 3, 1], [1, 1, -1]]
            sage: G3.number_of_blocks()
            2

        """
        return len(self._symbol)


    def determinant(self):
        """
        Returns the (`p`-part of the) determinant (square-class) of the
        Hessian matrix of the quadratic form (given by regarding the
        integral symmetric matrix which generated this genus symbol as
        the Gram matrix of `Q`) associated to this local genus symbol.

        OUTPUT:

        an integer

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            sage: G2.determinant()
            128

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G3
            Genus symbol at 3:     1^3 3^-1
            sage: G3.determinant()
            3
        """
        p = self._prime
        return prod([ p**(s[0]*s[1]) for s in self._symbol ])

    det = determinant

    def dimension(self):
        """
        Return the dimension of a quadratic form associated to this genus symbol.

        OUTPUT:

        an integer ``>= 0``

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            sage: G2.dimension()
            4

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 3
            sage: G3 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G3
            Genus symbol at 3:     1^3 3^-1
            sage: G3.dimension()
            4

        """
        return sum([ s[1] for s in self._symbol ])

    # TODO:
    # DELETE rank() IN FAVOR OF THE dimension() ... why?

    rank = dimension

    def direct_sum(self, other):
        r"""
        Return the local genus of the direct sum of two representatives.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
            sage: A = matrix.diagonal([1,2,3,4])
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [1^-2 2^1 4^1]_6
            sage: G2.direct_sum(G2)
            Genus symbol at 2:    [1^4 2^2, 4^2]_4
        """
        if self.prime() != other.prime():
            raise ValueError("the local genus symbols must be over the same prime")
        sym1 = self.symbol_tuple_list()
        sym2 = other.symbol_tuple_list()
        m = max(sym1[-1][0],sym2[-1][0])
        sym1 = dict([[s[0],s] for s in sym1])
        sym2 = dict([[s[0],s] for s in sym2])

        symbol = []
        for k in range(m+1):
            if self.prime()==2:
                b = [k, 0, 1, 0, 0]
            else:
                b = [k, 0, 1]
            for sym in [sym1, sym2]:
                try:
                    s = sym[k]
                    b[1] += s[1]
                    b[2] *= s[2]
                    if self.prime()==2:
                        b[2] = b[2] % 8
                        if s[3] == 1:
                            b[3] = s[3]
                        b[4] = (b[4] + s[4]) % 8
                except KeyError:
                    pass
            if b[1] != 0:
                symbol.append(b)
        return Genus_Symbol_p_adic_ring(self.prime(), symbol, check=True)

    def excess(self):
        """
        Returns the p-excess of the quadratic form whose Hessian
        matrix is the symmetric matrix A.  When p = 2 the p-excess is
        called the oddity.

        WARNING/NOTE: This normalization seems non-standard, and we
        should review this entire class to make sure that we have our
        doubling conventions straight throughout!

        REFERENCE:

        [Co1999]_ Conway and Sloane Book, 3rd edition, pp 370-371.

        OUTPUT:

        an integer

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: AC = diagonal_matrix(ZZ, [1,3,-3])
            sage: p=2; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            1
            sage: p=3; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=5; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=7; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=11; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0

            sage: AC = 2 * diagonal_matrix(ZZ, [1,3,-3])
            sage: p=2; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            1
            sage: p=3; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=5; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=7; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0
            sage: p=11; Genus_Symbol_p_adic_ring(p, p_adic_symbol(AC, p, 2)).excess()
            0

            sage: A = 2*diagonal_matrix(ZZ, [1,2,3,4])
            sage: p=2; Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)).excess()
            2
            sage: p=3; Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)).excess()
            6
            sage: p=5; Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)).excess()
            0
            sage: p=7; Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)).excess()
            0
            sage: p=11; Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)).excess()
            0

        """
        p = self._prime
        if self._prime == 2:
           k = 0
           for s in self._symbol:
               if s[0]%2 == 1 and s[2] in (3,5):
                   k += 1
           return Integer(sum([ s[4] for s in self._symbol ]) + 4*k).mod(8)
        else:
           k = 0
           for s in self._symbol:
               if s[0]%2 == 1 and s[2] == -1:
                   k += 1
           return Integer(sum([ s[1]*(p**s[0]-1) for s in self._symbol ]) + 4*k).mod(8)

    def scale(self):
        r"""
        Return the scale of this local genus.

        OUTPUT:

        an integer
        """
        return self.prime()**self._symbol[0][0]

    def level(self):
        r"""
        Return the maximal scale of a jordan component.
        """
        return self.prime()**self._symbol[-1][0]

    def trains(self):
        """
        Compute the indices for each of the trains in this local genus
        symbol if it is associated to the prime p=2 (and raise an
        error for all other primes).

        OUTPUT:

        a list of integers >= 0

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            sage: G2.trains()
            [[0, 1, 2]]

        """
        ## Check that p = 2
        if self._prime != 2:
            raise TypeError("trains() only makes sense when the prime of the p_adic_Genus_Symbol is p=2")
        symbol = self._symbol
        return canonical_2_adic_trains(symbol)


    def compartments(self):
        """
        Compute the indices for each of the compartments in this local genus
        symbol if it is associated to the prime p=2 (and raise an
        error for all other primes).

        OUTPUT:

        a list of integers >= 0

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import p_adic_symbol
            sage: from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: p = 2
            sage: G2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(A, p, 2)); G2
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            sage: G2.compartments()
            [[0, 1, 2]]

        """
        ## Check that p = 2
        if self._prime != 2:
            raise TypeError("compartments() only makes sense when the prime of the p_adic_Genus_Symbol is p=2")
        symbol = self._symbol
        return canonical_2_adic_compartments(symbol)

class GenusSymbol_global_ring(object):
    """
    This represents a collection of local genus symbols (at primes)
    and signature information which represent the genus of a
    non-degenerate integral lattice.

    INPUT:

    - ``signature_pair`` -- a tuple of two non-negative integers

    - ``local_symbols`` -- a list of :class:`Genus_Symbol_p_adic_ring`` instances

    - ``representative`` -- (default: ``None``) integer symmetric matrix
      the gram matrix of a representative of this genus

    - ``check`` -- (default: ``True``) a boolean; checks the input

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, LocalGenusSymbol
        sage: A = matrix.diagonal(ZZ, [2,4,6,8])
        sage: local_symbols = [LocalGenusSymbol(A, p) for p in (2*A.det()).prime_divisors()]
        sage: G = GenusSymbol_global_ring((4,0),local_symbols, representative=A);G
        Genus of
        [2 0 0 0]
        [0 4 0 0]
        [0 0 6 0]
        [0 0 0 8]
        Genus symbol at 2:    [2^-2 4^1 8^1]_6
        Genus symbol at 3:     1^3 3^-1
    """

    def __init__(self, signature_pair, local_symbols, representative=None, check=True):
        """
        Initialize a global genus symbol from a non-degenerate
        integral gram matrix (and possibly information about its
        largest elementary divisors).


        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus

            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: G = Genus(A)
            sage: G == loads(dumps(G))
            True
        """
        if check:
            if not all([type(sym)==Genus_Symbol_p_adic_ring for sym in local_symbols]):
                raise TypeError("local symbols must be a list of local genus symbols")
            n = signature_pair[0] + signature_pair[1]
            if not all([sym.dimension()==n for sym in local_symbols]):
                raise TypeError("all local symbols must be of the same dimension")
            if representative is not None:
                representative.is_symmetric()

        self._representative = representative
        self._signature = signature_pair
        self._local_symbols = local_symbols



    def __repr__(self):
        r"""
        Return a string representing the global genus symbol.

        OUTPUT:

        a string

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: A = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: GS = Genus(A)
            sage: GS
            Genus of
            [2 0 0 0]
            [0 4 0 0]
            [0 0 6 0]
            [0 0 0 8]
            Genus symbol at 2:    [2^-2 4^1 8^1]_6
            Genus symbol at 3:     1^3 3^-1

            sage: A2 = Matrix(ZZ,2,2,[2,-1,-1,2])
            sage: Genus(A2)
            Genus of
            [ 2 -1]
            [-1  2]
            Genus symbol at 2:    1^-2
            Genus symbol at 3:     1^-1 3^-1

        """
        local_symbols = "Signature: %s"%(self._signature,)
        for s in self._local_symbols:
            local_symbols += "\n" + s.__repr__()
        if self._representative is not None:
            return "Genus of\n%s\n%s" % (self._representative,local_symbols)
        else:
            return "Genus \n%s" %local_symbols


    def _latex_(self):
        """
        The Latex representation of this lattice.

        EXAMPLES::

            sage: D4=QuadraticForm(Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2]))
            sage: G=D4.global_genus_symbol()
            sage: G._latex_()
            '\\mbox{Genus of}\\\\\\left(\\begin{array}{rrrr}\n2 & 0 & 0 & -1 \\\\\n0 & 2 & 0 & -1 \\\\\n0 & 0 & 2 & -1 \\\\\n-1 & -1 & -1 & 2\n\\end{array}\\right)\\\\\\\\\\mbox{Genus symbol at } 2\\mbox{: }1^{-2}  :2^{-2} '
        """
        local_symbols = ""
        for s in self._local_symbols:
            local_symbols += "\\\\" + s._latex_()
        return "\\mbox{Genus of}\\\\%s\\\\%s" % (self._representative._latex_(),local_symbols)



    def __eq__(self, other):
        """
        Determines if two global genus symbols are equal (not just equivalent!).

        INPUT:

        a :class:`GenusSymbol_global_ring` object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus

            sage: A1 = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: GS1 = Genus(A1)
            sage: A2 = DiagonalQuadraticForm(ZZ, [1,2,3,5]).Hessian_matrix()
            sage: GS2 = Genus(A2)

            sage: GS1 == GS2
            False

            sage: GS2 == GS1
            False

            sage: GS1 == GS1
            True

            sage: GS2 == GS2
            True

        TESTS::

            sage: D4=QuadraticForm(Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2]))
            sage: G=D4.global_genus_symbol()
            sage: sage.quadratic_forms.genera.genus.is_GlobalGenus(G)
            True
            sage: G==deepcopy(G)
            True
            sage: sage.quadratic_forms.genera.genus.is_GlobalGenus(G)
            True
        """
        if self is other:
            return True
        t = len(self._local_symbols)
        if t != len(other._local_symbols):
            return False
        for i in range(t):
            if self._local_symbols[i] != other._local_symbols[i]:
                return False
        return True



    def __ne__(self, other):
        """
        Determine if two global genus symbols are unequal (not just inequivalent!).

        INPUT:

        a ``GenusSymbol_global_ring`` object

        OUTPUT:

        boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus

            sage: A1 = DiagonalQuadraticForm(ZZ, [1,2,3,4]).Hessian_matrix()
            sage: GS1 = Genus(A1)
            sage: A2 = DiagonalQuadraticForm(ZZ, [1,2,3,5]).Hessian_matrix()
            sage: GS2 = Genus(A2)

            sage: GS1 != GS2
            True

            sage: GS2 != GS1
            True

            sage: GS1 != GS1
            False

            sage: GS2 != GS2
            False

        """
        return not self == other

    def is_even(self):
        r"""
        Return if this genus is even.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: G = Genus(Matrix(ZZ,2,[2,1,1,2]))
            sage: G.is_even()
            True
        """
        return self._local_symbols[0].is_even()

    def signature_pair_of_matrix(self):
        """
        Return the signature pair (p, n) of the (non-degenerate)
        global genus symbol, where p is the number of positive
        eigenvalues and n is the number of negative eigenvalues.

        OUTPUT:

        a pair of integers `(p, n)` each `>= 0`

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus

            sage: A = DiagonalQuadraticForm(ZZ, [1,-2,3,4,8,-11]).Hessian_matrix()
            sage: GS = Genus(A)
            sage: GS.signature_pair_of_matrix()
            (4, 2)

        """
        return self._signature

    signature_pair = signature_pair_of_matrix

    def automorphous_numbers(self):
        r"""
        """
        aut_numbers = []
        for sym in self.local_symbols():
            aut_numbers.append((sym.prime(),sym.automorphous_numbers()))
        return aut_numbers

    def spinor_kernel(self):
        r"""
        """
        from sage.quadratic_forms.genera.spinor_genus import AdelicSquareClasses
        syms = self.local_symbols()
        primes = tuple([sym.prime() for sym in syms])
        grp = AdelicSquareClasses(primes)
        kernel_gens = []
        # -1 adic contribution
        sig = self.signature_pair_of_matrix()
        if sig[0]*sig[1] > 1:
            kernel_gens.append(grp.delta(-1, p=-1))
        for sym in syms:
            for A in sym.automorphous_numbers():
                kernel_gens.append(grp.delta(A, p=sym.prime()))
        return grp, grp.subgroup(kernel_gens)


    def determinant(self):
        """
        Return the determinant of this genus, where the determinant
        is the Hessian determinant of the quadratic form whose Gram
        matrix is the Gram matrix giving rise to this global genus
        symbol.

        OUTPUT:

        an integer

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: A = DiagonalQuadraticForm(ZZ, [1,-2,3,4]).Hessian_matrix()
            sage: GS = Genus(A)
            sage: GS.determinant()
            -384

        """
        r, s = self.signature_pair_of_matrix()
        return (-1)**s*prod([ G.determinant() for G in self._local_symbols ])

    det = determinant

    def dimension(self):
        r"""
        Return the dimension of this genus.
        """
        p, n = self.signature_pair_of_matrix()
        return p + n

    rank = dimension

    def discriminant_form(self):
        r"""
        Return the discriminant form associated to this genus.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: A = matrix.diagonal(ZZ, [2,-4,6,8])
            sage: GS = Genus(A)
            sage: GS.discriminant_form()
            Finite quadratic module over Integer Ring with invariants (2, 2, 4, 24)
            Gram matrix of the quadratic form with values in Q/2Z:
            [ 1/2    0    0    0]
            [   0  3/2    0    0]
            [   0    0  7/4    0]
            [   0    0    0 7/24]
            sage: A = matrix.diagonal(ZZ, [1,-4,6,8])
            sage: GS = Genus(A)
            sage: GS.discriminant_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 24)
            Gram matrix of the quadratic form with values in Q/Z:
            [ 1/2    0    0]
            [   0  3/4    0]
            [   0    0 7/24]
        """
        from sage.modules.torsion_quadratic_module import TorsionQuadraticForm
        qL = []
        for gs in self._local_symbols:
            p = gs._prime
            for block in gs.symbol_tuple_list():
                qL.append(_gram_from_jordan_block(p, block, True))
        q = matrix.block_diagonal(qL)
        return TorsionQuadraticForm(q)


    def direct_sum(self, other):
        r"""
        Return the genus of the direct sum of ``self`` and ``other``.

        The direct sum is defined as the direct sum of representatives.

        EXAMPLES::

            sage: G = IntegralLattice("A4").twist(3).genus()
            sage: G.direct_sum(G)
            Genus of
            None
            Genus symbol at 2:    1^8
            Genus symbol at 3:     3^8
            Genus symbol at 5:     1^6 5^2

        """
        p1, n1 = self.signature_pair()
        p2, n2 = other.signature_pair()
        signature_pair = (p1 + p2, n1 + n2)

        primes = [s.prime() for s in self.local_symbols()]
        primes += [s.prime() for s in other.local_symbols() if not s.prime() in primes]

        local_symbols = []
        for p in primes:
            sym_p = self.local_symbols(p=p).direct_sum(other.local_symbols(p=p))
            local_symbols.append(sym_p)
        return GenusSymbol_global_ring(signature_pair, local_symbols)

    def local_symbols(self, p=None):
        r"""
        Return the local symbols.

        INPUT:

        - ``p`` - a prime number; if given return only the symbol at `p`.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import Genus
            sage: A = matrix.diagonal(ZZ, [2,-4,6,8])
            sage: GS = Genus(A)
            sage: GS.local_symbols()
            [Genus symbol at 2:    [2^-2 4^1 8^1]_4,
             Genus symbol at 3:     1^-3 3^-1]

        """
        if p is None:
            return deepcopy(self._local_symbols)
        p = ZZ(p)
        for sym in self._local_symbols:
            if p == sym.prime():
                return deepcopy(sym)
        assert p!=2
        sym_p = [[0, self.rank(), self.det().kronecker(p)]]
        return Genus_Symbol_p_adic_ring(p, sym_p)

    def rational_representative(self):
        r"""

        """
        from sage.quadratic_forms.all import QuadraticForm
        sminus = self.signature_pair_of_matrix()[1]
        det = self.determinant()
        m = self.rank()
        P = []
        for sym in self._local_symbols:
            p = sym._prime
            # it is important to use the definition of cassels here!
            if QuadraticForm(QQ,2*sym.gram_matrix()).hasse_invariant(p) == -1:
                P.append(p)
        return rational_qf_from_invariants(m, det, P, sminus)

    def representative(self):
        r"""
        Return a representative in this genus.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.genus import all_genera_by_det, Genus
            sage: g = all_genera_by_det([1,3],24)[0]

        A representative is not know.
        Let us trigger its computation::

            sage: g == Genus(g.representative())
            True
        """
        if self._representative is None:
            self._compute_representative()
        return self._representative

    def level(self):
        r"""
        """
        return prod(sym.level() for sym in self.local_symbols())

    def scale(self):
        r"""
        Return the scale of this genus.

        OUTPUT:

        an integer

        EXAMPLES::

            sage: Genus(matrix.
        """
        return prod([s.scale() for s in self.local_symbols()])

    def _compute_representative(self, LLL=True):
        r"""
        Return a representative of this genus.


        """
        from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
        even = self.is_even()
        q = self.rational_representative()
        n = q.nrows()
        # the associated quadratic form xGx.T/2 should be integral
        L = IntegralLattice(4*q).maximal_overlattice().gram_matrix()
        p = 2
        sym2 = self.local_symbols()[0]
        if not self.is_even():
            # the quadratic form of xGx.T/2 must be integral
            # for things to work
            # solve this by multiplying the basis by 2
            L = local_modification(L, 4*sym2.gram_matrix(), p)
            L = (L/4).change_ring(ZZ)
        else:
            L = local_modification(L, sym2.gram_matrix(), p)
        for sym in self._local_symbols[1:]:
            p = sym.prime()
            L = local_modification(L, sym.gram_matrix(), p)
        # confirm the computation
        L = L.change_ring(ZZ)
        if LLL:
            sig = self.signature_pair_of_matrix()
            if sig[0]*sig[1] != 0:
                from sage.env import SAGE_EXTCODE
                from sage.interfaces.gp import Gp
                from cypari2.pari_instance import Pari
                pari = Pari()
                gp = Gp()
                m = pari(L)
                gp.read(SAGE_EXTCODE + "/pari/simon/qfsolve.gp")
                m = gp.eval('qflllgram_indefgoon(%s)'%m)
                # convert the output string to sage
                L = pari(m).sage()[0]
            elif sig[1] != 0:
                U = -(-L).LLL_gram()
                L = U.T * L * U
            else:
                U = L.LLL_gram()
                L = U.T * L * U
        assert Genus(L) == self
        self._representative = L

    def representatives(self, algorithm="default"):
        r"""
        Return a list of representatives for the classes in this genus
        """
        n = self.dimension()
        representatives = []
        if algorithm == "default":
            if n == 2 and ZZ.prod(self.signature_pair_of_matrix())!=0:
                # binary indefinite
                algorithm = "sage"
            else:
                algorithm = "magma"
        if algorithm == "magma":
            from sage.interfaces.magma import Magma
            magma = Magma()
            magma.set_server_and_command(command="/LOCAL/magma/magma")
            if prod(self.signature_pair_of_matrix()) != 0:
                if n <= 2:
                    raise NotImplementedError()
                K = magma.RationalsAsNumberField()
                gram = magma.Matrix(K, n, self.representative().list())
                L = gram.NumberFieldLatticeWithGram()
                representatives = L.GenusRepresentatives()
                representatives = [r.GramMatrix().ChangeRing(magma.Rationals()).sage() for r in representatives]
                return representatives
            else:
                e = 1
                if self.signature_pair_of_matrix()[1] != 0:
                    e = -1
                K = magma.Rationals()
                gram = magma.Matrix(K, n, (e*self.representative()).list())
                L = gram.LatticeWithGram()
                representatives = L.GenusRepresentatives()
                representatives = [e*r.GramMatrix().sage() for r in representatives]
                return representatives

        elif algorithm == "sage":
            if n > 2:
                raise NotImplementedError()
            if n == 1:
                return self.representative()
            if n == 2:
                d = - 4 * self.determinant()
                from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
                for q in BinaryQF_reduced_representatives(d, proper=False):
                    if q[1] % 2 == 0:  # we want integrality of the gram matrix
                        m = matrix(ZZ, 2, [q[0], q[1]//2, q[1]//2, q[2]])
                        if Genus(m) == self:
                            representatives.append(m)
                # recompute using magma
                if ZZ.prod(self.signature_pair_of_matrix()) == 0:
                    assert len(representatives)==len(self.representatives(algorithm="magma"))
        else:
            raise ValueError("unknown algorithm")
        return representatives




def rational_qf_from_invariants(m, det, P, sminus):
    r"""
    Return the gram matrix of a rational quadratic form with given invariants.

    INPUT:

    - ``m`` -- integer; the rank

    - ``det`` -- rational; the determinant

    - ``P`` -- a list of primes where cassels hasse invariant is negative

    - ``sminus`` -- integer; negative part of the signature

    OUTPUT:

    - a symmetric matrix; the gram matrix


    """
    from sage.arith.misc import hilbert_symbol
    from sage.rings.infinity import Infinity
    P = [ZZ(p) for p in P]
    m = ZZ(m)
    d = ZZ(det).squarefree_part()
    sminus = ZZ(sminus)
    if d.sign() != (-1)**sminus:
        raise ValueError("Invariants do not define a rational quadratic form")
    if m == 1 and len(P) != 0:
        raise ValueError("Invariants do not define a rational quadratic form")
    if m == 2:
        for p in P:
            if QQ(-d).is_padic_square(p):
                raise ValueError("Invariants do not define a rational quadratic form")
    if sminus % 4 in (2, 3) and len(P) % 2 == 0:
        raise ValueError("")
    D = []
    while m >= 2:
        if m >= 4:
            if sminus > 0:
                a = ZZ(-1)
            else:
                a = ZZ(1)
        elif m == 3:
            Pprime = [p for p in P if hilbert_symbol(-1, -d, p)==1]
            Pprime += [p for p in (2*d).prime_divisors() if hilbert_symbol(-1, -d, p)==-1 and p not in P]
            if sminus > 0:
                a = ZZ(-1)
            else:
                a = ZZ(1)
            for p in Pprime:
                if d.valuation(p) % 2 == 0:
                    a *= p
            assert all([(a*d).valuation(p)%2==1 for p in Pprime])
        elif m == 2:
            S = P
            if sminus == 2:
                S += [-1]
            a = QQ.hilbert_conductor_inverse(S,-d)
        P = [p for p in P if hilbert_symbol(a, -d, p) == 1] + [p for p in ZZ(2*a*d).prime_divisors() if hilbert_symbol(a, -d, p)==-1 and p not in P]
        sminus = max(0, sminus-1)
        m = m - 1
        d = a*d
        D.append(ZZ(a).squarefree_part())
    d = ZZ(d).squarefree_part()
    D.append(d.squarefree_part())
    return matrix.diagonal(D)

def local_modification(M, Lp, p, check=True):
    r"""
    Return a local modification of `M` that matches `Lp` at `p`.

    INPUT:

    - ``M`` -- the gram matrix of a `ZZp`-maximal lattice

    - ``Lp`` -- gram matrix of a lattice with Lp `\QQ_p` rationally equivalent
      to `M`

    - ``p`` -- a prime number

    OUTPUT:

    a gram matrix `M'` such that `M` and `M'` are locally equivalent at all
    primes except `p` where `M'` is locally equivalent to `Lp`

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import local_modification, Genus
        sage: L = IntegralLattice("A3").twist(15)
        sage: M = L.maximal_overlattice().gram_matrix()
        sage: L = L.gram_matrix()
        sage: for p in prime_divisors(L.det()):
        ....:     M = local_modification(M, L, p)
        sage: Genus(M) == Genus(L)
        True
        sage: from sage.quadratic_forms.genera.genus import local_modification
        sage: L = IntegralLattice("D4").twist(3*4)
        sage: M = L.maximal_overlattice()
        sage: local_modification(M.gram_matrix(), L.gram_matrix(), 2)
        [16  8  8 16]
        [ 8  8  4 12]
        [ 8  4 24  0]
        [16 12  0 24]
    """
    from sage.quadratic_forms.genera.normal_form import p_adic_normal_form
    from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
    from sage.rings.finite_rings.finite_field_constructor import GF
    from sage.matrix.special import random_matrix
    M = matrix(M)
    Lp = matrix(Lp) # target lattice
    d = Lp.inverse().denominator()
    n = M.rank()
    scale = d.valuation(p)
    d = p**scale

    M = IntegralLattice(M)
    Lp_max = IntegralLattice(Lp).maximal_overlattice(p=p)

    # invert the gerstein operations
    _, U = p_adic_normal_form(Lp_max.gram_matrix(), p, precision=scale+3)
    B = (~Lp_max.basis_matrix()).change_ring(ZZ)*~U.change_ring(ZZ)

    _, UM = p_adic_normal_form(M.gram_matrix(), p, precision=scale+3)
    B = B * UM.change_ring(ZZ) * M.basis_matrix()

    # the local modification
    S = M.sublattice(((M.span(B) & M) + d * M).gens())
    S = S.gram_matrix()
    # confirm result
    if check:
        s1 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(S, p, scale))
        s2 = Genus_Symbol_p_adic_ring(p, p_adic_symbol(Lp, p, scale))
        assert s1 == s2, "oops"
    return S

def _gram_from_jordan_block(p, block, discr_form=False):
    r"""
    Return the gram matrix of this jordan block.

    This is a helper for :meth:`discriminant_form` and :meth:`gram_matrix`.
    No input checks.

    INPUT:

    - ``p`` -- a prime number

    - ``block`` -- a list of 3 integers or 5 integers if `p` is `2`

    - ``discr_form`` -- bool (default: ``False``); if ``True`` invert the scales
      to obtain a gram matrix for the discriminant form instead.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.genus import _gram_from_jordan_block
        sage: block = [1, 3, 1]
        sage: _gram_from_jordan_block(5, block)
        [5 0 0]
        [0 5 0]
        [0 0 5]
        sage: block = [1, 4, 7, 1, 2]
        sage: _gram_from_jordan_block(2, block)
        [0 2 0 0]
        [2 0 0 0]
        [0 0 2 0]
        [0 0 0 2]

    For the discriminant form we obtain::

        sage: block = [1, 3, 1]
        sage: _gram_from_jordan_block(5, block, True)
        [4/5   0   0]
        [  0 2/5   0]
        [  0   0 2/5]
        sage: block = [1, 4, 7, 1, 2]
        sage: _gram_from_jordan_block(2, block, True)
        [  0 1/2   0   0]
        [1/2   0   0   0]
        [  0   0 1/2   0]
        [  0   0   0 1/2]
    """
    from sage.quadratic_forms.genera.normal_form import _min_nonsquare
    scale = block[0]
    rk = block[1]
    det = block[2]
    if p == 2:
        o = ZZ(block[3])
        t = ZZ(block[4])
        U = matrix(QQ,2,[0,1,1,0])
        V = matrix(QQ,2,[2,1,1,2])
        W = matrix(QQ,1,[1])
        if o == 0:
            if det in [1, 7]:
                qL = (rk // 2) * [U]
            else:
                qL = (rk//2 - 1)*[U] + [V]
        if o == 1:
            if rk % 2 == 1:
                qL = max(0, (rk - 3) // 2) * [U]
                if t*det % 8 in [3, 5]:
                    qL += [V]
                elif rk >= 3:
                    qL += [U]
                qL += [t * W]
            else:
                if det in [3, 5]:
                    det = -1
                else:
                    det = 1
                qL = max(0, (rk - 4) // 2) * [U]
                if (det , t) == (1, 0):
                    qL += [U, 1 * W, 7 * W]
                if (det , t) == (1, 2):
                    qL += [U, 1 * W, 1 * W]
                if (det , t) == (1, 4):
                    qL += [V, 1 * W, 3 * W]
                if (det , t) == (1, 6):
                    qL += [U, 7 * W, 7 * W]
                if (det , t) == (-1, 0):
                    qL += [V, 1 * W, 7 * W]
                if (det , t) == (-1, 2):
                    qL += [U, 3 * W, 7 * W]
                if (det , t) == (-1, 4):
                    qL += [U, 1 * W, 3 * W]
                if (det , t) == (-1, 6):
                    qL += [U, 1 * W, 5 * W]
                # if the rank is 2 there is a U too much
                if rk == 2:
                    qL = qL[-2:]
        q = matrix.block_diagonal(qL)
        if discr_form:
            q = q / 2**scale
        else:
            q = q * 2**scale
    if p != 2 and discr_form:
        q = matrix.identity(QQ, rk)
        d = 2**(rk % 2)
        if Integer(d).kronecker(p) != det:
            u = ZZ(_min_nonsquare(p))
            q[0,0] = u
        q = q * (2 / p**scale)
    if p != 2 and not discr_form:
        q = matrix.identity(QQ, rk)
        if det != 1:
            u = ZZ(_min_nonsquare(p))
            q[0,0] = u
        q = q * p**scale
    return q
