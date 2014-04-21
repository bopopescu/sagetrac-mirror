from sage.misc.cachefunc import cached_method
from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2


class ZinbielOperad(CombinatorialFreeModule):
    r"""
    The Zinbiel operad

    EXAMPLES::

        sage: Z = ZinbielOperad(QQ)
        sage: Z.an_element()
        B[word: ]
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = ZinbielOperad(QQ); A
            The Zinbiel operad over Rational Field
            sage: TestSuite(A).run()

            sage: W = Words()
            sage: A.composition(A(W("abc")), A(W("def")), "b")
            B[word: adcef] + B[word: adecf] + B[word: adefc]
            sage: A.composition(A("abc"), A("def"), "b")
            B[word: adcef] + B[word: adecf] + B[word: adefc]
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category=OperadsWithBasis(R))

    def _repr_(self):
        """
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: ZinbielOperad(QQ)     # indirect doctest
            The Zinbiel operad over Rational Field
        """
        return "The Zinbiel operad over {}".format(self.base_ring())

    def species(self):
        """
        Returns the species of non-empty lists.

        This is the species associated with the Zinbiel operad.

        EXAMPLES::

            sage: f = ZinbielOperad(QQ).species()
            sage: f.generating_series().coefficients(5)
            [0, 1, 1, 1, 1]
        """
        from sage.combinat.species.library import LinearOrderSpecies
        return LinearOrderSpecies().restricted(min=1)

    def _coerce_end(self, st):
        """
        Allows for the shortcut ``A(<string>)``

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: A("abc")
            B[word: abc]
        """
        if isinstance(st, str):
            return self._from_key(st)
        raise TypeError

    def _from_key(self, k):
        """
        Returns an element from a word.

        EXAMPLES::

            sage: ZinbielOperad(QQ)._from_key("abc")
            B[word: abc]
        """
        return self._element_constructor(self.basis().keys()(k))

    @cached_method
    def one_basis(self, letter=1):
        """
        Returns the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def operad_generators(self):
        """
        EXAMPLES::

            sage: ZinbielOperad(QQ).operad_generators()
            Finite family {'zinbiel_product': B[word: 12]}
        """
        from sage.sets.family import Family
        return Family(dict([("zinbiel_product",
                             self._from_key([1, 2]))]))

    def composition_on_basis_list(self, x, y, i):
        r"""
        Returns the composition of two words `x o_i y` as a list of
        words.

        The composition index `i` must be a label of `x`.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis_list(Words("abc"), Words("de"),"a")
            [word: dbce, word: dbec, word: debc]
        """
        if not(i in x):
            raise ValueError("the composition index is not present")
        elif x[0] == i:
            return map(lambda u: y[:1] + u, ShuffleProduct_w1w2(x[1:], y[1:]))
        return map(lambda u: x[:1] + u,
                   self.composition_on_basis_list(x[1:], y, i))

    def composition_on_basis(self, x, y, i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words(["a","b","c"]), Words(["d","e"]),"a")
            B[word: dbce] + B[word: dbec] + B[word: debc]
        """
        if not(i in x):
            raise ValueError("the composition index is not present")
        return self.sum_of_monomials(t for t in
                                     self.composition_on_basis_list(x, y, i))
