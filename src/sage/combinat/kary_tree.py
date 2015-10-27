# -*- coding: utf-8 -*-
"""
K-ary Trees

This module deals with k-ary trees as mathematical (in particular immutable)
objects.

AUTHORS:

- Adrien Boussicault, Bérénice Delcroix-Oger (2015): initial implementation.

"""
#*****************************************************************************
#       Copyright (C) 2015 
#           Adrien Boussicault <boussica@labri.fr>, 
#           Bérénice Delcroix-Oger <berenice.delcroix@math.univ-toulouse.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.list_clone import ClonableArray
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.combinat.ordered_tree import LabelledOrderedTrees
from sage.rings.integer import Integer
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.combinat.combinatorial_map import combinatorial_map
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.misc.cachefunc import cached_method


class KaryTree(AbstractClonableTree, ClonableArray):
    """
    K-ary trees.

    K-ary trees here mean ordered (a.k.a. plane) finite k-ary
    trees, where "ordered" means that the children of each node are
    ordered.

    K-ary trees contain nodes and leaves, where each node has k
    children while each leaf has no children. The number of leaves
    of a k-ary tree always equals k-1 times the number of nodes, plus `1`.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length k of k-ary trees or convertible objects. This corresponds
      to the standard recursive definition of a k-ary tree as either a
      leaf or a k-tuple of k-ary trees. Syntactic sugar allows leaving out
      all but the outermost calls of the ``KaryTree()`` constructor, so
      that, e. g., ``KaryTree([KaryTree(None),KaryTree(None)])`` can
      be shortened to ``KaryTree([None,None])``. The arity is given by the 
      size of the list given in parameter.

    - `̀̀̀`arity`` -- ``None`` (default) or a positive integer. This corresponds 
      to the arity of the tree. If ``None`` is given then the constructor will 
      try to deduce the arity from the size of ``children``.

    - ``check`` -- (default: ``True``) whether check for k-arity should be
      performed or not.

    EXAMPLES::

        sage: KaryTree(None)
        .
        sage: KaryTree(None,3)
        .
        sage: KaryTree([None, None, None])
        [., ., .]
        sage: KaryTree([None, None, None], 3)
        [., ., .]
        sage: KaryTree([None, [None, None], None])
        [., [., .], .]
        sage: KaryTree([[None, None], None])
        [[., .], .]
        sage: KaryTree("[[., ., .], ., .]")
        [[., ., .], ., .]
        sage: KaryTree([None, KaryTree([None, None])])
        [., [., .]]
        sage: KaryTree([KaryTree([None])])
        [[.]]

        sage: KaryTree([[None, None], None, None])
        Traceback (most recent call last):
        ...
        ValueError: this is not a k-ary tree

    TESTS::

        #sage: t1 = BinaryTree([[None, [[],[[], None]]],[[],[]]])
        #sage: t2 = BinaryTree([[[],[]],[]])
        #sage: with t1.clone() as t1c:
        #....:     t1c[1,1,1] = t2
        #sage: t1 == t1c
        #False
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        TODO DOC
        Ensure that binary trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`BinaryTree`.

        TESTS::

#            sage: from sage.combinat.binary_tree import BinaryTrees_all
#            sage: issubclass(BinaryTrees_all().element_class, BinaryTree)
#            True
#            sage: t0 = BinaryTree([[],[[], None]])
#            sage: t0.parent()
#            Binary trees
#            sage: type(t0)
#            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>
#
#            sage: t1 = BinaryTrees()([[],[[], None]])
#            sage: t1.parent() is t0.parent()
#            True
#            sage: type(t1) is type(t0)
#            True
#
#            sage: t1 = BinaryTrees(4)([[],[[], None]])
#            sage: t1.parent() is t0.parent()
#            True
#            sage: type(t1) is type(t0)
#            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        TODO DOC
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

#            sage: BinaryTree._auto_parent
#            Binary trees
#            sage: BinaryTree([None, None]).parent()
#            Binary trees
        """
        return KaryTrees_all()

    def __init__(self, parent, children = None, arity = None, check = True):
        """
        TESTS::

            sage: KaryTree([None, None]).parent()
            K-ary trees
            sage: KaryTree("[., [., [., .]]]")
            [., [., [., .]]]
            sage: KaryTree("[., [., ., .], .]")
            [., [., ., .], .]
            sage: KaryTree("[.,.,.]", 4)
            Traceback (most recent call last):
            ...
            ValueError: wrong arity
#            sage: all(KaryTree(repr(t)) == t for n in range(6) for t in KaryTrees(k,n))
#            True
        """
        if (isinstance(children, str)):  # if the input is the repr of a binary tree
            children = children.replace(".","None")
            from ast import literal_eval
            children = literal_eval(children)
        if children is None:
            children = []
        elif (children == [] or children == ()
              or isinstance(children, (Integer, int))):
            children = [None, None]
        if (children.__class__ is self.__class__ and
            children.parent() == parent):
            children = list(children)
        else:
            children = [self.__class__(parent, x, len( children ) ) for x in children]
        if arity is None:
            self._arity = len( children )
        else:
            self._arity = arity
        ClonableArray.__init__(self, parent, children, check=check)

    def arity(self):
        return len(self)

    def check(self):
        """
        Check that ``self`` is a k-ary tree.

        EXAMPLES::

            sage: KaryTree([[None, None], [None, None]])     # indirect doctest
            [[., .], [., .]]
            sage: BinaryTree(None) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a k-ary tree
        """
        if self and len(self) != self._arity :
            raise ValueError("this is not a k-ary tree")

    def _repr_(self):
        """
        TESTS::

            sage: t1 = BinaryTree([[], None]); t1  # indirect doctest
            [[., .], .]
            sage: BinaryTree([[None, t1], None])   # indirect doctest
            [[., [[., .], .]], .]
        """
        if not self:
            return "."
        else:
            return super(KaryTree, self)._repr_()


    def is_empty(self):
        """
        Return whether ``self`` is empty.

        The notion of emptiness employed here is the one which defines
        a binary tree to be empty if its root is a leaf. There is
        precisely one empty k-ary tree.

        EXAMPLES::

            sage: KaryTree().is_empty()
            True
            sage: KaryTree([]).is_empty()
            False
            sage: KaryTree([[], None]).is_empty()
            False
        """
        return not self


    def make_leaf(self):
        """
        Modify ``self`` so that it becomes a leaf (i. e., an empty tree).

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::

            :meth:`make_node <sage.combinat.binary_tree.BinaryTree.make_node>`

        EXAMPLES::

            sage: t = KaryTree([None, None])
            sage: t.make_leaf()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ....:     t1.make_leaf()
            sage: t, t1
            ([., .], .)
        """
        self._require_mutable()
        self.__init__(self.parent(), None)

    def comb(self, side=0):
        r"""
        Return the comb of a k-ary trees : there are k combs in a k-ary
        tree (one for every direction). A comb is defined as the list of  
        children of the nodes on branches whose direction does not change.

        INPUT:
        - ``side`` -- set to 0 to obtain the leftmost comb, to i+1 to
        obtain the leftmost comb in the set of combs which are on the
        right side of the ith comb and to k-1 to obtain the rightmost comb
        of a k-ary tree.

        OUTPUT:

        A list of $k$ tuples of k-ary trees.

        EXAMPLES::

            sage: T = KaryTree([[None,None, None],[[None,None, [None, None, None]], None, None],None])
            sage: T.comb(0)
            [[None, ., .]]
            sage: T.comb(1)
            [[[., ., [., ., .]], None, .]]
            sage: T.comb(2)
            []
            sage: T = KaryTree([[[[None]]]])
            sage: T.comb(0)
            [[None], [None], [None]]


        """
        if self.is_empty():
            return []
        d=self.arity()
        if not side < d :
            raise ValueError("Value %d is a wrong side value : it must be strictly smaller than the arity %d of the tree"%(side,d)) 
        tree=self[side]
        res=[]
        fc=[]
        while not tree.is_empty():
            for i in range(d):
                if i==side:
                    fc.append(None)
                else:
                    fc.append(tree[i])
            res.append(fc)
            fc=[]
            tree=tree[side]
        return res

    def hook_number(self):
        r"""
        Return the number of hooks in a k-ary trees.

        The hook of a vertex v is the union of {v}, and all the  
        branches from {v} in which the direction does not change. 

        There is a unique way to partition the vertices in hooks.   
        The number of hooks in such the partition is the hook number 
        of the tree.

        We can obtain this partition recursively by extracting the root's
        hook and iterating the processus on each tree of the remaining
        forest.      

        EXAMPLES::
            sage: T = KaryTree(None)
            sage: T.hook_number()
            0     
            sage: T = KaryTree( [None,None,None] )
            sage: T.hook_number()
            1
            sage: T = KaryTree([[None, [None, None]], [[None, None], None]])
            sage: T.hook_number()
            3
            sage: T = KaryTree( [None,None,None] )
            sage: T.hook_number()
            1
        """
        if self.is_empty() or self==None:
            return 0
        s=1
        for i in range(self.arity()):
            for h in self.comb(i):
                if len(h)>0:
                    for el in h: 
                        if not(el==None) and not(el.is_empty()):
                            s+=el.hook_number()
        return s

#    def _ascii_art_( self ):
#        r"""
#        TESTS::
#
#            sage: ascii_art(BinaryTree())
#            <BLANKLINE>
#            sage: ascii_art(BinaryTree([]))
#            o
#            sage: for bt in BinaryTrees(3):
#            ....:     print ascii_art(bt)
#            o
#             \
#              o
#               \
#                o
#            o
#             \
#              o
#             /
#            o
#              o
#             / \
#            o   o
#              o
#             /
#            o
#             \
#              o
#                o
#               /
#              o
#             /
#            o
#            sage: ascii_art(BinaryTree([None,[]]))
#            o
#             \
#              o
#            sage: ascii_art(BinaryTree([None,[None,[]]]))
#            o
#             \
#              o
#               \
#                o
#            sage: ascii_art(BinaryTree([None,[[],None]]))
#            o
#             \
#              o
#             /
#            o
#            sage: ascii_art(BinaryTree([None,[[[],[]],[]]]))
#               o
#                \
#                _o_
#               /   \
#              o     o
#             / \
#            o   o
#            sage: ascii_art(BinaryTree([None,[[None,[[],[]]],None]]))
#            o
#             \
#              o
#             /
#            o
#             \
#              o
#             / \
#            o   o
#            sage: ascii_art(BinaryTree([[],None]))
#              o
#             /
#            o
#            sage: ascii_art(BinaryTree([[[[],None], None],None]))
#                  o
#                 /
#                o
#               /
#              o
#             /
#            o
#            sage: ascii_art(BinaryTree([[[],[]],None]))
#                o
#               /
#              o
#             / \
#            o   o
#            sage: ascii_art(BinaryTree([[[None,[]],[[[],None],None]], None]))
#                   o
#                  /
#              ___o___
#             /       \
#            o         o
#             \       /
#              o     o
#                   /
#                  o
#            sage: ascii_art(BinaryTree([[None,[[],[]]],None]))
#              o
#             /
#            o
#             \
#              o
#             / \
#            o   o
#            sage: ascii_art(BinaryTree([[],[]]))
#              o
#             / \
#            o   o
#            sage: ascii_art(BinaryTree([[],[[],None]]))
#              _o_
#             /   \
#            o     o
#                 /
#                o
#            sage: ascii_art(BinaryTree([[None,[]],[[[],None],None]]))
#              ___o___
#             /       \
#            o         o
#             \       /
#              o     o
#                   /
#                  o
#            sage: ascii_art(BinaryTree([[[],[]],[[],None]]))
#                __o__
#               /     \
#              o       o
#             / \     /
#            o   o   o
#            sage: ascii_art(BinaryTree([[[],[]],[[],[]]]))
#                __o__
#               /     \
#              o       o
#             / \     / \
#            o   o   o   o
#            sage: ascii_art(BinaryTree([[[[],[]],[[],[]]],[]]))
#                    ___o___
#                   /       \
#                __o__       o
#               /     \
#              o       o
#             / \     / \
#            o   o   o   o
#            sage: ascii_art(BinaryTree([[],[[[[],[]],[[],[]]],[]]]))
#              _____o______
#             /            \
#            o           ___o___
#                       /       \
#                    __o__       o
#                   /     \
#                  o       o
#                 / \     / \
#                o   o   o   o
#        """
#        node_to_str = lambda bt: str(bt.label()) if hasattr(bt, "label") else "o"
#
#        if self.is_empty():
#            from sage.typeset.ascii_art import empty_ascii_art
#            return empty_ascii_art
#
#        from sage.typeset.ascii_art import AsciiArt
#        if self[0].is_empty() and self[1].is_empty():
#            bt_repr = AsciiArt( [node_to_str(self)] )
#            bt_repr._root = 1
#            return bt_repr
#        if self[0].is_empty():
#            node = node_to_str(self)
#            rr_tree = self[1]._ascii_art_()
#            if rr_tree._root > 2:
#                f_line = " " ** Integer( rr_tree._root - 3 ) + node
#                s_line = " " ** Integer( len( node ) + rr_tree._root - 3 ) + "\\"
#                t_repr = AsciiArt( [f_line, s_line] ) * rr_tree
#                t_repr._root = rr_tree._root - 2
#            else:
#                f_line = node
#                s_line = " " + "\\"
#                t_line = " " ** Integer( len( node ) + 1 )
#                t_repr = AsciiArt( [f_line, s_line] ) * ( AsciiArt( [t_line] ) + rr_tree )
#                t_repr._root = rr_tree._root
#            t_repr._baseline = t_repr._h - 1
#            return t_repr
#        if self[1].is_empty():
#            node = node_to_str(self)
#            lr_tree = self[0]._ascii_art_()
#            f_line = " " ** Integer( lr_tree._root + 1 ) + node
#            s_line = " " ** Integer( lr_tree._root ) + "/"
#            t_repr = AsciiArt( [f_line, s_line] ) * lr_tree
#            t_repr._root = lr_tree._root + 2
#            t_repr._baseline = t_repr._h - 1
#            return t_repr
#        node = node_to_str(self)
#        lr_tree = self[0]._ascii_art_()
#        rr_tree = self[1]._ascii_art_()
#        nb_ = lr_tree._l - lr_tree._root + rr_tree._root - 1
#        nb_L = int( nb_ / 2 )
#        nb_R = nb_L + ( 1 if nb_ % 2 == 1 else 0 )
#        f_line = " " ** Integer( lr_tree._root + 1 ) + "_" ** Integer( nb_L ) + node
#        f_line += "_" ** Integer( nb_R )
#        s_line = " " ** Integer( lr_tree._root ) + "/" + " " ** Integer( len( node ) + rr_tree._root - 1 + ( lr_tree._l - lr_tree._root ) ) + "\\"
#        t_repr = AsciiArt( [f_line, s_line] ) * ( lr_tree + AsciiArt( [" " ** Integer( len( node ) + 2 )] ) + rr_tree )
#        t_repr._root = lr_tree._root + nb_L + 2
#        t_repr._baseline = t_repr._h - 1
#        return t_repr
#
#    def is_empty(self):
#        """
#        Return whether ``self`` is empty.
#
#        The notion of emptiness employed here is the one which defines
#        a binary tree to be empty if its root is a leaf. There is
#        precisely one empty binary tree.
#
#        EXAMPLES::
#
#            sage: BinaryTree().is_empty()
#            True
#            sage: BinaryTree([]).is_empty()
#            False
#            sage: BinaryTree([[], None]).is_empty()
#            False
#        """
#        return not self
#
#    def graph(self, with_leaves=True):
#        """
#        Convert ``self`` to a digraph. By default, this graph contains
#        both nodes and leaves, hence is never empty. To obtain a graph
#        which contains only the nodes, the ``with_leaves`` optional
#        keyword variable has to be set to ``False``.
#
#        INPUT:
#
#        - ``with_leaves`` -- (default: ``True``) a Boolean, determining
#          whether the resulting graph will be formed from the leaves
#          and the nodes of ``self`` (if ``True``), or only from the
#          nodes of ``self`` (if ``False``)
#
#        EXAMPLES::
#
#            sage: t1 = BinaryTree([[], None])
#            sage: t1.graph()
#            Digraph on 5 vertices
#            sage: t1.graph(with_leaves=False)
#            Digraph on 2 vertices
#
#            sage: t1 = BinaryTree([[], [[], None]])
#            sage: t1.graph()
#            Digraph on 9 vertices
#            sage: t1.graph().edges()
#            [(0, 1, None), (0, 4, None), (1, 2, None), (1, 3, None), (4, 5, None), (4, 8, None), (5, 6, None), (5, 7, None)]
#            sage: t1.graph(with_leaves=False)
#            Digraph on 4 vertices
#            sage: t1.graph(with_leaves=False).edges()
#            [(0, 1, None), (0, 2, None), (2, 3, None)]
#
#            sage: t1 = BinaryTree()
#            sage: t1.graph()
#            Digraph on 1 vertex
#            sage: t1.graph(with_leaves=False)
#            Digraph on 0 vertices
#
#            sage: BinaryTree([]).graph()
#            Digraph on 3 vertices
#            sage: BinaryTree([]).graph(with_leaves=False)
#            Digraph on 1 vertex
#
#            sage: t1 = BinaryTree([[], [[], []]])
#            sage: t1.graph(with_leaves=False)
#            Digraph on 5 vertices
#            sage: t1.graph(with_leaves=False).edges()
#            [(0, 1, None), (0, 2, None), (2, 3, None), (2, 4, None)]
#        """
#        from sage.graphs.graph import DiGraph
#
#        if with_leaves:   # We want leaves and nodes.
#
#            # Special treatment for the case when self is empty.
#            # In this case, rec(self, 0) would give a false result.
#            if not self:
#                return DiGraph({0: []})
#
#            res = DiGraph()
#            # The edge set of res will be built up step by step using the
#            # following function:
#            def rec(tr, idx):
#                if not tr:  # tr is a leaf.
#                    return
#                else:  # tr is a node.
#                    nbl = 2 * tr[0].node_number() + 1
#                    res.add_edges([[idx, idx + 1], [idx, idx + 1 + nbl]])
#                    rec(tr[0], idx + 1)
#                    rec(tr[1], idx + nbl + 1)
#            rec(self, 0)
#            return res
#
#        else:   # We want only the nodes.
#
#            # Special treatment for the case when self has only 1 node.
#            # In this case, the general DiGraph construction would
#            # falsely yield an empty graph (since it adds nodes only
#            # implicitly by adding edges).
#            if self.node_number() == 1:
#                return DiGraph({0: []})
#
#            res = DiGraph()
#            # The edge set of res will be built up step by step using the
#            # following function:
#            def rec(tr, idx):
#                if not tr:  # tr is a leaf.
#                    return
#                else:  # tr is a node.
#                    nbl = tr[0].node_number()
#                    if nbl > 0:
#                        res.add_edge([idx, idx + 1])
#                        rec(tr[0], idx + 1)
#                    if tr[1].node_number() > 0:
#                        res.add_edge([idx, idx + nbl + 1])
#                        rec(tr[1], idx + nbl + 1)
#            rec(self, 0)
#            return res
#
#    def canonical_labelling(self, shift=1):
#        r"""
#        Return a labelled version of ``self``.
#
#        The canonical labelling of a binary tree is a certain labelling of the
#        nodes (not the leaves) of the tree.
#        The actual canonical labelling is currently unspecified. However, it
#        is guaranteed to have labels in `1...n` where `n` is the number of
#        nodes of the tree. Moreover, two (unlabelled) trees compare as equal if
#        and only if their canonical labelled trees compare as equal.
#
#        EXAMPLES::
#
#            sage: BinaryTree().canonical_labelling()
#            .
#            sage: BinaryTree([]).canonical_labelling()
#            1[., .]
#            sage: BinaryTree([[[], [[], None]], [[], []]]).canonical_labelling()
#            5[2[1[., .], 4[3[., .], .]], 7[6[., .], 8[., .]]]
#        """
#        LTR = self.parent().labelled_trees()
#        if self:
#            sz0 = self[0].node_number()
#            return LTR([self[0].canonical_labelling(shift),
#                        self[1].canonical_labelling(shift+1+sz0)],
#                       label=shift+sz0)
#        else:
#            return LTR(None)
#
#    def show(self, with_leaves=False):
#        """
#        Show the binary tree ``show``, with or without leaves depending
#        on the Boolean keyword variable ``with_leaves``.
#
#        .. WARNING::
#
#            Left and right children might get interchanged in
#            the actual picture. Moreover, for a labelled binary
#            tree, the labels shown in the picture are not (in
#            general) the ones given by the labelling!
#
#            Use :meth:`_latex_`, ``view``,
#            :meth:`_ascii_art_` or ``pretty_print`` for more
#            faithful representations of the data of the tree.
#
#        TESTS::
#
#            sage: t1 = BinaryTree([[], [[], None]])
#            sage: t1.show()
#        """
#        try:
#            self.graph(with_leaves=with_leaves).show(layout='tree', tree_root=0, tree_orientation="down")
#        except RuntimeError:
#            # This is for the border case BinaryTree().show().
#            self.graph(with_leaves=with_leaves).show()
#
#    def make_node(self, child_list = [None, None]):
#        """
#        Modify ``self`` so that it becomes a node with children ``child_list``.
#
#        INPUT:
#
#        - ``child_list`` -- a pair of binary trees (or objects convertible to)
#
#        .. NOTE:: ``self`` must be in a mutable state.
#
#        .. SEEALSO::
#
#            :meth:`make_leaf <sage.combinat.binary_tree.BinaryTree.make_leaf>`
#
#        EXAMPLES::
#
#            sage: t = BinaryTree()
#            sage: t.make_node([None, None])
#            Traceback (most recent call last):
#            ...
#            ValueError: object is immutable; please change a copy instead.
#            sage: with t.clone() as t1:
#            ....:     t1.make_node([None, None])
#            sage: t, t1
#            (., [., .])
#            sage: with t.clone() as t:
#            ....:     t.make_node([BinaryTree(), BinaryTree(), BinaryTree([])])
#            Traceback (most recent call last):
#            ...
#            ValueError: the list must have length 2
#            sage: with t1.clone() as t2:
#            ....:     t2.make_node([t1, t1])
#            sage: with t2.clone() as t3:
#            ....:     t3.make_node([t1, t2])
#            sage: t1, t2, t3
#            ([., .], [[., .], [., .]], [[., .], [[., .], [., .]]])
#        """
#        self._require_mutable()
#        child_lst = [self.__class__(self.parent(), x) for x in child_list]
#        if not(len(child_lst) == 2):
#            raise ValueError("the list must have length 2")
#        self.__init__(self.parent(), child_lst, check=False)
#
#    def make_leaf(self):
#        """
#        Modify ``self`` so that it becomes a leaf (i. e., an empty tree).
#
#        .. NOTE:: ``self`` must be in a mutable state.
#
#        .. SEEALSO::
#
#            :meth:`make_node <sage.combinat.binary_tree.BinaryTree.make_node>`
#
#        EXAMPLES::
#
#            sage: t = BinaryTree([None, None])
#            sage: t.make_leaf()
#            Traceback (most recent call last):
#            ...
#            ValueError: object is immutable; please change a copy instead.
#            sage: with t.clone() as t1:
#            ....:     t1.make_leaf()
#            sage: t, t1
#            ([., .], .)
#        """
#        self._require_mutable()
#        self.__init__(self.parent(), None)
#
#    def _to_ordered_tree(self, bijection="left", root=None):
#        r"""
#        Internal recursive method to obtain an ordered tree from a binary
#        tree.
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([[],[]])
#            sage: bt._to_ordered_tree()
#            [[], [[]]]
#            sage: bt._to_ordered_tree(bijection="right")
#            [[[]], []]
#            sage: bt._to_ordered_tree(bijection="none")
#            Traceback (most recent call last):
#            ...
#            ValueError: the bijection argument should be either left or right
#            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
#            sage: bt._to_ordered_tree()
#            [[], [[], []], [[], [[]]]]
#            sage: bt._to_ordered_tree(bijection="right")
#            [[[[]], [[]]], [[]], []]
#        """
#        close_root = False
#        if(root is None):
#            from sage.combinat.ordered_tree import OrderedTree
#            root = OrderedTree().clone()
#            close_root = True
#        if(self):
#            left, right = self[0],self[1]
#            if(bijection == "left"):
#                root = left._to_ordered_tree(bijection=bijection,root=root)
#                root.append(right._to_ordered_tree(bijection=bijection,root=None))
#            elif(bijection =="right"):
#                root.append(left._to_ordered_tree(bijection=bijection, root=None))
#                root = right._to_ordered_tree(bijection=bijection,root=root)
#            else:
#                raise ValueError("the bijection argument should be either left or right")
#        if(close_root):
#            root.set_immutable()
#        return root
#
#    @combinatorial_map(name="To ordered tree, left child = left brother")
#    def to_ordered_tree_left_branch(self):
#        r"""
#        Return an ordered tree of size `n+1` by the following recursive rule:
#
#        - if `x` is the left child of `y`, `x` becomes the left brother
#          of `y`
#        - if `x` is the right child of `y`, `x` becomes the last child
#          of `y`
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([[],[]])
#            sage: bt.to_ordered_tree_left_branch()
#            [[], [[]]]
#            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
#            sage: bt.to_ordered_tree_left_branch()
#            [[], [[], []], [[], [[]]]]
#        """
#        return self._to_ordered_tree()
#
#    @combinatorial_map(name="To ordered tree, right child = right brother")
#    def to_ordered_tree_right_branch(self):
#        r"""
#        Return an ordered tree of size `n+1` by the following recursive rule:
#
#        - if `x` is the right child of `y`, `x` becomes the right brother
#          of `y`
#        - if `x` is the left child of `y`, `x` becomes the first child
#          of `y`
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([[],[]])
#            sage: bt.to_ordered_tree_right_branch()
#            [[[]], []]
#            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
#            sage: bt.to_ordered_tree_right_branch()
#            [[[[]], [[]]], [[]], []]
#        """
#        return self._to_ordered_tree(bijection="right")
#
#    def _postfix_word(self, left_first = True, start = 1):
#        r"""
#        Internal recursive method to obtain a postfix canonical read of the
#        binary tree.
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([[],[]])
#            sage: bt._postfix_word()
#            [1, 3, 2]
#            sage: bt._postfix_word(left_first=False)
#            [3, 1, 2]
#            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
#            sage: bt._postfix_word()
#            [1, 3, 4, 2, 6, 8, 7, 5]
#            sage: bt._postfix_word(left_first=False)
#            [8, 6, 7, 3, 4, 1, 2, 5]
#        """
#        if not self:
#            return []
#        left = self[0]._postfix_word(left_first, start)
#        label = start + self[0].node_number()
#        right = self[1]._postfix_word(left_first, start = label +1)
#        if left_first:
#            left.extend(right)
#            left.append(label)
#            return left
#        else:
#            right.extend(left)
#            right.append(label)
#            return right
#
#    @combinatorial_map(name="To complete tree")
#    def as_ordered_tree(self, with_leaves=True):
#        r"""
#        Return the same tree seen as an ordered tree. By default, leaves
#        are transformed into actual nodes, but this can be avoided by
#        setting the optional variable ``with_leaves`` to ``False``.
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([]); bt
#            [., .]
#            sage: bt.as_ordered_tree()
#            [[], []]
#            sage: bt.as_ordered_tree(with_leaves = False)
#            []
#            sage: bt = bt.canonical_labelling(); bt
#            1[., .]
#            sage: bt.as_ordered_tree()
#            1[None[], None[]]
#            sage: bt.as_ordered_tree(with_leaves=False)
#            1[]
#        """
#        if with_leaves:
#            children = [child.as_ordered_tree(with_leaves) for child in self]
#        else:
#            if not self:
#                raise ValueError("The empty binary tree cannot be made into an ordered tree with with_leaves = False")
#            children = [child.as_ordered_tree(with_leaves) for child in self if not child.is_empty()]
#        if self in LabelledBinaryTrees():
#            from sage.combinat.ordered_tree import LabelledOrderedTree
#            return LabelledOrderedTree(children, label = self.label())
#        else:
#            from sage.combinat.ordered_tree import OrderedTree
#            return OrderedTree(children)
#
#    @combinatorial_map(name="To graph")
#    def to_undirected_graph(self, with_leaves=False):
#        r"""
#        Return the undirected graph obtained from the tree nodes and edges.
#
#        Leaves are ignored by default, but one can set ``with_leaves`` to
#        ``True`` to obtain the graph of the complete tree.
#
#        INPUT:
#
#        - ``with_leaves`` -- (default: ``False``) a Boolean, determining
#          whether the resulting graph will be formed from the leaves
#          and the nodes of ``self`` (if ``True``), or only from the
#          nodes of ``self`` (if ``False``)
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([])
#            sage: bt.to_undirected_graph()
#            Graph on 1 vertex
#            sage: bt.to_undirected_graph(with_leaves=True)
#            Graph on 3 vertices
#
#            sage: bt = BinaryTree()
#            sage: bt.to_undirected_graph()
#            Graph on 0 vertices
#            sage: bt.to_undirected_graph(with_leaves=True)
#            Graph on 1 vertex
#
#        If the tree is labelled, we use its labelling to label the graph.
#        Otherwise, we use the graph canonical labelling which means that
#        two different trees can have the same graph.
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([[],[None,[]]])
#            sage: bt.canonical_labelling().to_undirected_graph() == bt.to_undirected_graph()
#            False
#            sage: BinaryTree([[],[]]).to_undirected_graph() == BinaryTree([[[],None],None]).to_undirected_graph()
#            True
#        """
#        if (not with_leaves) and (not self):
#            # this case needs extra care :(
#            from sage.graphs.graph import Graph
#            return Graph([])
#        return self.as_ordered_tree(with_leaves).to_undirected_graph()
#
#    @combinatorial_map(name="To poset")
#    def to_poset(self, with_leaves=False, root_to_leaf=False):
#        r"""
#        Return the poset obtained by interpreting the tree as a Hasse
#        diagram.
#
#        The default orientation is from leaves to root but you can
#        pass ``root_to_leaf=True`` to obtain the inverse orientation.
#
#        Leaves are ignored by default, but one can set ``with_leaves`` to
#        ``True`` to obtain the poset of the complete tree.
#
#        INPUT:
#
#        - ``with_leaves`` -- (default: ``False``) a Boolean, determining
#          whether the resulting poset will be formed from the leaves
#          and the nodes of ``self`` (if ``True``), or only from the
#          nodes of ``self`` (if ``False``)
#        - ``root_to_leaf`` -- (default: ``False``) a Boolean,
#          determining whether the poset orientation should be from root
#          to leaves (if ``True``) or from leaves to root (if ``False``).
#
#        EXAMPLES::
#
#            sage: bt = BinaryTree([])
#            sage: bt.to_poset()
#            Finite poset containing 1 elements
#            sage: bt.to_poset(with_leaves=True)
#            Finite poset containing 3 elements
#            sage: P1 = bt.to_poset(with_leaves=True)
#            sage: len(P1.maximal_elements())
#            1
#            sage: len(P1.minimal_elements())
#            2
#            sage: bt = BinaryTree([])
#            sage: P2 = bt.to_poset(with_leaves=True,root_to_leaf=True)
#            sage: len(P2.maximal_elements())
#            2
#            sage: len(P2.minimal_elements())
#            1
#
#        If the tree is labelled, we use its labelling to label the poset.
#        Otherwise, we use the poset canonical labelling::
#
#            sage: bt = BinaryTree([[],[None,[]]]).canonical_labelling()
#            sage: bt
#            2[1[., .], 3[., 4[., .]]]
#            sage: bt.to_poset().cover_relations()
#            [[4, 3], [3, 2], [1, 2]]
#
#        Let us check that the empty binary tree is correctly handled::
#
#            sage: bt = BinaryTree()
#            sage: bt.to_poset()
#            Finite poset containing 0 elements
#            sage: bt.to_poset(with_leaves=True)
#            Finite poset containing 1 elements
#        """
#        if (not with_leaves) and (not self):
#            # this case needs extra care :(
#            from sage.combinat.posets.posets import Poset
#            return Poset({})
#        return self.as_ordered_tree(with_leaves).to_poset(root_to_leaf)
#
#    @combinatorial_map(order = 2, name="Left-right symmetry")
#    def left_right_symmetry(self):
#        r"""
#        Return the left-right symmetrized tree of ``self``.
#
#        EXAMPLES::
#
#            sage: BinaryTree().left_right_symmetry()
#            .
#            sage: BinaryTree([]).left_right_symmetry()
#            [., .]
#            sage: BinaryTree([[],None]).left_right_symmetry()
#            [., [., .]]
#            sage: BinaryTree([[None, []],None]).left_right_symmetry()
#            [., [[., .], .]]
#        """
#        if not self:
#            return BinaryTree()
#        tree = [self[1].left_right_symmetry(),self[0].left_right_symmetry()]
#        if(not self in LabelledBinaryTrees()):
#            return BinaryTree(tree)
#        return LabelledBinaryTree(tree, label = self.label())
#
#    @combinatorial_map(order=2, name="Left border symmetry")
#    def left_border_symmetry(self):
#        r"""
#        Return the tree where a symmetry has been applied recursively on
#        all left borders. If a tree is made of three trees `[T_1, T_2,
#        T_3]` on its left border, it becomes `[T_3', T_2', T_1']` where
#        same symmetry has been applied to `T_1, T_2, T_3`.
#
#        EXAMPLES::
#
#            sage: BinaryTree().left_border_symmetry()
#            .
#            sage: BinaryTree([]).left_border_symmetry()
#            [., .]
#            sage: BinaryTree([[None,[]],None]).left_border_symmetry()
#            [[., .], [., .]]
#            sage: BinaryTree([[None,[None,[]]],None]).left_border_symmetry()
#            [[., .], [., [., .]]]
#            sage: bt = BinaryTree([[None,[None,[]]],None]).canonical_labelling()
#            sage: bt
#            4[1[., 2[., 3[., .]]], .]
#            sage: bt.left_border_symmetry()
#            1[4[., .], 2[., 3[., .]]]
#        """
#        if not self:
#            return BinaryTree()
#        border = []
#        labelled = self in LabelledBinaryTrees()
#        labels = []
#        t = self
#        while(t):
#            border.append(t[1].left_border_symmetry())
#            if labelled: labels.append(t.label())
#            t = t[0]
#        tree = BinaryTree()
#        for r in border:
#            if labelled:
#                tree = LabelledBinaryTree([tree,r],label=labels.pop(0))
#            else:
#                tree = BinaryTree([tree,r])
#        return tree
#
#    def in_order_traversal_iter(self):
#        """
#        The depth-first infix-order traversal iterator for the binary
#        tree ``self``.
#
#        This method iters each vertex (node and leaf alike) of the given
#        binary tree following the depth-first infix order traversal
#        algorithm.
#
#        The *depth-first infix order traversal algorithm* iterates
#        through a binary tree as follows::
#
#            iterate through the left subtree (by the depth-first infix
#                order traversal algorithm);
#            yield the root;
#            iterate through the right subtree (by the depth-first infix
#                order traversal algorithm).
#
#        For example on the following binary tree `T`, where we denote
#        leaves by `a, b, c, \ldots` and nodes by `1, 2, 3, \ldots`::
#
#            |     ____3____          |
#            |    /         \         |
#            |   1          __7__     |
#            |  / \        /     \    |
#            | a   2      _5_     8   |
#            |    / \    /   \   / \  |
#            |   b   c  4     6 h   i |
#            |         / \   / \      |
#            |        d   e f   g     |
#
#        the depth-first infix-order traversal algorithm iterates through
#        the vertices of `T` in the following order:
#        `a,1,b,2,c,3,d,4,e,5,f,6,g,7,h,8,i`.
#
#        See :meth:`in_order_traversal` for a version of this algorithm
#        which not only iterates through, but actually does something at
#        the vertices of tree.
#
#        TESTS::
#
#            sage: b = BinaryTree([[],[[],[]]]); ascii_art([b])
#            [   _o_     ]
#            [  /   \    ]
#            [ o     o   ]
#            [      / \  ]
#            [     o   o ]
#            sage: ascii_art(list(b.in_order_traversal_iter()))
#            [ , o, ,   _o_    , , o, ,   o  , , o,  ]
#            [         /   \             / \         ]
#            [        o     o           o   o        ]
#            [             / \                       ]
#            [            o   o                      ]
#            sage: ascii_art(filter(lambda node: node.label() is not None,
#            ....:     b.canonical_labelling().in_order_traversal_iter()))
#            [ 1,   _2_    , 3,   4  , 5 ]
#            [     /   \         / \     ]
#            [    1     4       3   5    ]
#            [         / \               ]
#            [        3   5              ]
#
#            sage: list(BinaryTree(None).in_order_traversal_iter())
#            [.]
#        """
#        if self.is_empty():
#            yield self
#            return
#        # TODO:: PYTHON 3
#        # yield from self[0].in_order_traversal_iter()
#        for left_subtree in self[0].in_order_traversal_iter():
#            yield left_subtree
#        yield self
#        # TODO:: PYTHON 3
#        # yield from self[1].in_order_traversal_iter()
#        for right_subtree in self[1].in_order_traversal_iter():
#            yield right_subtree
#
#    def in_order_traversal(self, node_action=None, leaf_action=None):
#        r"""
#        Explore the binary tree ``self`` using the depth-first infix-order
#        traversal algorithm, executing the ``node_action`` function
#        whenever traversing a node and executing the ``leaf_action``
#        function whenever traversing a leaf.
#
#        In more detail, what this method does to a tree `T` is the
#        following::
#
#            if the root of `T` is a node:
#                apply in_order_traversal to the left subtree of `T`
#                    (with the same node_action and leaf_action);
#                apply node_action to the root of `T`;
#                apply in_order_traversal to the right subtree of `T`
#                    (with the same node_action and leaf_action);
#            else:
#                apply leaf_action to the root of `T`.
#
#        For example on the following binary tree `T`, where we denote
#        leaves by `a, b, c, \ldots` and nodes by `1, 2, 3, \ldots`::
#
#            |     ____3____          |
#            |    /         \         |
#            |   1          __7__     |
#            |  / \        /     \    |
#            | a   2      _5_     8   |
#            |    / \    /   \   / \  |
#            |   b   c  4     6 h   i |
#            |         / \   / \      |
#            |        d   e f   g     |
#
#        this method first applies ``leaf_action`` to `a`, then applies
#        ``node_action`` to `1`, then ``leaf_action`` to `b`, then
#        ``node_action`` to `2`, etc., with the vertices being traversed
#        in the order `a,1,b,2,c,3,d,4,e,5,f,6,g,7,h,8,i`.
#
#        See :meth:`in_order_traversal_iter` for a version of this
#        algorithm which only iterates through the vertices rather than
#        applying any function to them.
#
#        INPUT:
#
#        - ``node_action`` -- (optional) a function which takes a node in input
#          and does something during the exploration
#        - ``leaf_action`` -- (optional) a function which takes a leaf in input
#          and does something during the exploration
#
#        TESTS::
#
#            sage: nb_leaf = 0
#            sage: def l_action(_):
#            ....:    global nb_leaf
#            ....:    nb_leaf += 1
#            sage: nb_node = 0
#            sage: def n_action(_):
#            ....:    global nb_node
#            ....:    nb_node += 1
#
#            sage: BinaryTree().in_order_traversal(n_action, l_action)
#            sage: nb_leaf, nb_node
#            (1, 0)
#
#            sage: nb_leaf, nb_node = 0, 0
#            sage: b = BinaryTree([[],[[],[]]]); b
#            [[., .], [[., .], [., .]]]
#            sage: b.in_order_traversal(n_action, l_action)
#            sage: nb_leaf, nb_node
#            (6, 5)
#            sage: nb_leaf, nb_node = 0, 0
#            sage: b = b.canonical_labelling()
#            sage: b.in_order_traversal(n_action, l_action)
#            sage: nb_leaf, nb_node
#            (6, 5)
#            sage: l = []
#            sage: b.in_order_traversal(lambda node: l.append( node.label() ))
#            sage: l
#            [1, 2, 3, 4, 5]
#
#            sage: leaf = 'a'
#            sage: l = []
#            sage: def l_action(_):
#            ....:    global leaf, l
#            ....:    l.append(leaf)
#            ....:    leaf = chr( ord(leaf)+1 )
#            sage: n_action = lambda node: l.append( node.label() )
#            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
#            ....:     canonical_labelling()
#            sage: b.in_order_traversal(n_action, l_action)
#            sage: l
#            ['a', 1, 'b', 2, 'c', 3, 'd', 4, 'e', 5, 'f', 6, 'g', 7, 'h', 8,
#             'i']
#        """
#        if leaf_action is None:
#            leaf_action = lambda x: None
#        if node_action is None:
#            node_action = lambda x: None
#
#        for node in self.in_order_traversal_iter():
#            if node.is_empty():
#                leaf_action(node)
#            else:
#                node_action(node)
#
#    @combinatorial_map(name="Over operation on Binary Trees")
#    def over(self, bt):
#        r"""
#        Return ``self`` over ``bt``, where "over" is the ``over``
#        (`/`) operation.
#
#        If `T` and `T'` are two binary trees, then `T` over `T'`
#        (written `T / T'`) is defined as the tree obtained by grafting
#        `T'` on the rightmost leaf of `T`. More precisely, `T / T'` is
#        defined by identifying the root of the `T'` with the rightmost
#        leaf of `T`. See section 4.5 of [HNT05]_.
#
#        If `T` is empty, then `T / T' = T'`.
#
#        The definition of this "over" operation goes back to
#        Loday-Ronco [LodRon0102066]_ (Definition 2.2), but it is
#        denoted by `\backslash` and called the "under" operation there.
#        In fact, trees in sage have their root at the top, contrary to
#        the trees in [LodRon0102066]_ which are growing upwards. For
#        this reason, the names of the over and under operations are
#        swapped, in order to keep a graphical meaning.
#        (Our notation follows that of section 4.5 of [HNT05]_.)
#
#        .. SEEALSO::
#
#            :meth:`under`
#
#        EXAMPLES:
#
#        Showing only the nodes of a binary tree, here is an
#        example for the over operation::
#
#            |   o       __o__       _o_         |
#            |  / \  /  /     \  =  /   \        |
#            | o   o   o       o   o     o       |
#            |          \     /           \      |
#            |           o   o           __o__   |
#            |                          /     \  |
#            |                         o       o |
#            |                          \     /  |
#            |                           o   o   |
#
#        A Sage example::
#
#            sage: b1 = BinaryTree([[],[[],[]]])
#            sage: b2 = BinaryTree([[None, []],[]])
#            sage: ascii_art((b1, b2, b1/b2))
#            (   _o_    ,   _o_  ,   _o_           )
#            (  /   \      /   \    /   \          )
#            ( o     o    o     o  o     o_        )
#            (      / \    \            /  \       )
#            (     o   o    o          o    o      )
#            (                               \     )
#            (                               _o_   )
#            (                              /   \  )
#            (                             o     o )
#            (                              \      )
#            (                               o     )
#
#        TESTS::
#
#            sage: b1 = BinaryTree([[],[]]); ascii_art([b1])
#            [   o   ]
#            [  / \  ]
#            [ o   o ]
#            sage: b2 = BinaryTree([[None,[]],[[],None]]); ascii_art([b2])
#            [   __o__   ]
#            [  /     \  ]
#            [ o       o ]
#            [  \     /  ]
#            [   o   o   ]
#            sage: ascii_art([b1.over(b2)])
#            [   _o_         ]
#            [  /   \        ]
#            [ o     o       ]
#            [        \      ]
#            [       __o__   ]
#            [      /     \  ]
#            [     o       o ]
#            [      \     /  ]
#            [       o   o   ]
#
#        The same in the labelled case::
#
#            sage: b1 = b1.canonical_labelling()
#            sage: b2 = b2.canonical_labelling()
#            sage: ascii_art([b1.over(b2)])
#            [   _2_         ]
#            [  /   \        ]
#            [ 1     3       ]
#            [        \      ]
#            [       __3__   ]
#            [      /     \  ]
#            [     1       5 ]
#            [      \     /  ]
#            [       2   4   ]
#        """
#        B = self.parent()._element_constructor_
#        if self.is_empty():
#            return bt
#        if hasattr(self, "label"):
#            lab = self.label()
#            return B([self[0], self[1].over(bt)], lab)
#        else:
#            return B([self[0], self[1].over(bt)])
#
#    __div__ = over
#
#    @combinatorial_map(name="Under operation on Binary Trees")
#    def under(self, bt):
#        r"""
#        Return ``self`` under ``bt``, where "under" is the ``under``
#        (`\backslash`) operation.
#
#        If `T` and `T'` are two binary trees, then `T` under `T'`
#        (written `T \backslash T'`) is defined as the tree obtained
#        by grafting `T` on the leftmost leaf of `T'`. More precisely,
#        `T \backslash T'` is defined by identifying the root of `T`
#        with the leftmost leaf of `T'`.
#
#        If `T'` is empty, then `T \backslash T' = T`.
#
#        The definition of this "under" operation goes back to
#        Loday-Ronco [LodRon0102066]_ (Definition 2.2), but it is
#        denoted by `/` and called the "over" operation there. In fact,
#        trees in sage have their root at the top, contrary to the trees
#        in [LodRon0102066]_ which are growing upwards. For this reason,
#        the names of the over and under operations are swapped, in
#        order to keep a graphical meaning.
#        (Our notation follows that of section 4.5 of [HNT05]_.)
#
#        .. SEEALSO::
#
#            :meth:`over`
#
#        EXAMPLES:
#
#        Showing only the nodes of a binary tree, here is an
#        example for the under operation::
#
#            sage: b1 = BinaryTree([[],[]])
#            sage: b2 = BinaryTree([None,[]])
#            sage: ascii_art((b1, b2, b1 \ b2))
#            (   o  , o  ,     _o_   )
#            (  / \    \      /   \  )
#            ( o   o    o    o     o )
#            (              / \      )
#            (             o   o     )
#
#        TESTS::
#
#            sage: b1 = BinaryTree([[],[[None,[]],None]]); ascii_art([b1])
#            [   _o_   ]
#            [  /   \  ]
#            [ o     o ]
#            [      /  ]
#            [     o   ]
#            [      \  ]
#            [       o ]
#            sage: b2 = BinaryTree([[],[None,[]]]); ascii_art([b2])
#            [   o     ]
#            [  / \    ]
#            [ o   o   ]
#            [      \  ]
#            [       o ]
#            sage: ascii_art([b1.under(b2)])
#            [        o_     ]
#            [       /  \    ]
#            [      o    o   ]
#            [     /      \  ]
#            [   _o_       o ]
#            [  /   \        ]
#            [ o     o       ]
#            [      /        ]
#            [     o         ]
#            [      \        ]
#            [       o       ]
#
#        The same in the labelled case::
#
#            sage: b1 = b1.canonical_labelling()
#            sage: b2 = b2.canonical_labelling()
#            sage: ascii_art([b1.under(b2)])
#            [        2_     ]
#            [       /  \    ]
#            [      1    3   ]
#            [     /      \  ]
#            [   _2_       4 ]
#            [  /   \        ]
#            [ 1     5       ]
#            [      /        ]
#            [     3         ]
#            [      \        ]
#            [       4       ]
#        """
#        B = self.parent()._element_constructor_
#        if bt.is_empty():
#            return self
#        lab = None
#        if hasattr(bt, "label"):
#            lab = bt.label()
#            return B([self.under(bt[0]), bt[1]], lab)
#        else:
#            return B([self.under(bt[0]), bt[1]])
#
#    _backslash_ = under
#
#    def is_full(self):
#        r"""
#        Return ``True`` if ``self`` is full, else return ``False``.
#
#        A full binary tree is a tree in which every node either has two
#        child nodes or has two child leaves.
#
#        This is also known as *proper binary tree* or *2-tree* or *strictly
#        binary tree*.
#
#        For example::
#
#            |       __o__   |
#            |      /     \  |
#            |     o       o |
#            |    / \        |
#            |   o   o       |
#            |  /     \      |
#            | o       o     |
#
#        is not full but the next one is::
#
#            |         ___o___   |
#            |        /       \  |
#            |     __o__       o |
#            |    /     \        |
#            |   o       o       |
#            |  / \     / \      |
#            | o   o   o   o     |
#
#        EXAMPLES::
#
#            sage: BinaryTree([[[[],None],[None,[]]], []]).is_full()
#            False
#            sage: BinaryTree([[[[],[]],[[],[]]], []]).is_full()
#            True
#            sage: ascii_art(filter(lambda bt: bt.is_full(), BinaryTrees(5)))
#            [   _o_    ,     _o_   ]
#            [  /   \        /   \  ]
#            [ o     o      o     o ]
#            [      / \    / \      ]
#            [     o   o  o   o     ]
#        """
#        if self.is_empty():
#            return True
#        if self[0].is_empty() != self[1].is_empty():
#            return False
#        return self[0].is_full() and self[1].is_full()
#
#    def is_perfect(self):
#        r"""
#        Return ``True`` if ``self`` is perfect, else return ``False``.
#
#        A perfect binary tree is a full tree in which all leaves are at the
#        same depth.
#
#        For example::
#
#            |         ___o___   |
#            |        /       \  |
#            |     __o__       o |
#            |    /     \        |
#            |   o       o       |
#            |  / \     / \      |
#            | o   o   o   o     |
#
#        is not perfect but the next one is::
#
#            |     __o__     |
#            |    /     \    |
#            |   o       o   |
#            |  / \     / \  |
#            | o   o   o   o |
#
#        EXAMPLES::
#
#            sage: lst = lambda i: filter(lambda bt: bt.is_perfect(), BinaryTrees(i))
#            sage: for i in range(10): ascii_art(lst(i)) # long time
#            [  ]
#            [ o ]
#            [  ]
#            [   o   ]
#            [  / \  ]
#            [ o   o ]
#            [  ]
#            [  ]
#            [  ]
#            [     __o__     ]
#            [    /     \    ]
#            [   o       o   ]
#            [  / \     / \  ]
#            [ o   o   o   o ]
#            [  ]
#            [  ]
#        """
#        return 2 ** self.depth() - 1 == self.node_number()
#
#    def is_complete(self):
#        r"""
#        Return ``True`` if ``self`` is complete, else return ``False``.
#
#        In a nutshell, a complete binary tree is a perfect binary tree
#        except possibly in the last level, with all nodes in the last
#        level "flush to the left".
#
#        In more detail:
#        A complete binary tree (also called binary heap) is a binary tree in
#        which every level, except possibly the last one (the deepest), is
#        completely filled. At depth `n`, all nodes must be as far left as
#        possible.
#
#        For example::
#
#            |         ___o___   |
#            |        /       \  |
#            |     __o__       o |
#            |    /     \        |
#            |   o       o       |
#            |  / \     / \      |
#            | o   o   o   o     |
#
#        is not complete but the following ones are::
#
#            |     __o__          _o_            ___o___     |
#            |    /     \        /   \          /       \    |
#            |   o       o      o     o      __o__       o   |
#            |  / \     / \    / \          /     \     / \  |
#            | o   o   o   o, o   o    ,   o       o   o   o |
#            |                            / \     /          |
#            |                           o   o   o           |
#
#        EXAMPLES::
#
#            sage: lst = lambda i: filter(lambda bt: bt.is_complete(), BinaryTrees(i))
#            sage: for i in range(9): ascii_art(lst(i)) # long time
#            [  ]
#            [ o ]
#            [   o ]
#            [  /  ]
#            [ o   ]
#            [   o   ]
#            [  / \  ]
#            [ o   o ]
#            [     o   ]
#            [    / \  ]
#            [   o   o ]
#            [  /      ]
#            [ o       ]
#            [     _o_   ]
#            [    /   \  ]
#            [   o     o ]
#            [  / \      ]
#            [ o   o     ]
#            [     __o__   ]
#            [    /     \  ]
#            [   o       o ]
#            [  / \     /  ]
#            [ o   o   o   ]
#            [     __o__     ]
#            [    /     \    ]
#            [   o       o   ]
#            [  / \     / \  ]
#            [ o   o   o   o ]
#            [       __o__     ]
#            [      /     \    ]
#            [     o       o   ]
#            [    / \     / \  ]
#            [   o   o   o   o ]
#            [  /              ]
#            [ o               ]
#        """
#        if self.is_empty():
#            return True
#        # self := L ^ R
#        dL = self[0].depth()
#        dR = self[1].depth()
#        # if L is perfect
#        if self[0].is_perfect():
#            # if the depth of R == depth of L then R must be complete
#            if dL == dR:
#                return self[1].is_complete()
#            # else R is perfect with depth equals depth of L - 1
#            elif dL == dR + 1:
#                return self[1].is_perfect()
#            return False
#        # L is not perfect then R is perfect and the depth of L = the depth of
#        # R + 1
#        return self[0].is_complete() and self[1].is_perfect() and dL == dR + 1


# Abstract class to serve as a Factory no instance are created.
class KaryTrees(UniqueRepresentation, Parent):
    """
    Factory for k-ary trees.

    INPUT:

    - ``arity`` -- (optional) an integer
    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all k-ary trees (of the given ``arity`` and ``size`` if 
      specified )

    EXAMPLES::

        sage: KaryTrees()
        k-ary trees

        sage: KaryTrees(3)
        3-ary trees

        sage: KaryTrees(3, 2)
        3-ary trees of size 2

    .. NOTE:: this is a factory class whose constructor returns instances of
              subclasses.

    .. NOTE:: the fact that KaryTrees is a class instead of a simple callable
              is an implementation detail. It could be changed in the future
              and one should not rely on it.
    """
    @staticmethod
    def __classcall_private__(cls, arity=None, size=None):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import KaryTrees_all, KaryTrees_size
            sage: isinstance(KaryTrees(3, 2), KaryTrees)
            True
            sage: isinstance(KaryTrees(3), KaryTrees)
            True
            sage: isinstance(KaryTrees(), KaryTrees)
            True
            sage: KaryTrees(3, 2) is KaryTrees_size(3,2)
            True
            sage: KaryTrees(2, 5).cardinality()
            42
            sage: KaryTrees(3) is KaryTrees_all()
            True
        """
        if size is None and arity is None:
            return KaryTrees_all()
        else:
            if not (isinstance(arity, (Integer, int)) and arity >= 0):
                raise ValueError("arity must be a positive integer")
            if size is None:
                return KaryTrees_arity(Integer(arity))
            else:
                if not (isinstance(size, (Integer, int)) and size >= 0):
                    raise ValueError("size must be a nonnegative integer")
                return KaryTrees_size(Integer(arity), Integer(size))

    @cached_method
    def leaf(self):
        """
        Return a leaf tree with ``self`` as parent.

        EXAMPLES::

            sage: KaryTrees(3).leaf()
            .

        TEST::

            sage: (KaryTrees(3).leaf() is
            ....:  sage.combinat.kary_tree.KaryTrees_all(3).leaf())
            True
        """
        return self(None)

    def arity(self):
        r"""
        Returns the arity.

        EXAMPLES::

            sage: KaryTrees(3).arity()
            3
        """
        return self._arity




from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer
from sage.structure.list_clone import ClonableArray


class IntegerPair(ClonableArray):
    def __init__(self, parent, v, check = True):
        ClonableArray.__init__(self, parent, v, check=check)
    def check(self):
        if not (not self or len(self) == 2):
            raise ValueError("this is not a pair o integer")
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)
    @lazy_class_attribute
    def _auto_parent(cls):
        return IntegerPairs()


class IntegerPairs(UniqueRepresentation, Parent):
    def __init__(self):
        Parent.__init__(self, category = InfiniteEnumeratedSets())

    def _repr_(self):
        return "Pair of non negative integers"

    def __contains__(self, elt):
        return Integer(elt[0]) >= Integer(0) and Integer(elt[0]) >= Integer(0) 

    def __iter__(self):
        n = 0
        while True:
            for i in range(n+1):
                yield self._element_constructor_([i, n-i])
            n += 1

    def __call__(self, elt):
        if elt in self:
            return self._element_constructor_( elt )
        else:
            raise ValueError("Value %s is not a pair of non negative integer."%(elt))

    def an_element(self):
        return self._element_constructor_([Integer(42),Integer(3)])

    def next(self, v):
        if v[1] == 0:
            return self._element_constructor_( [0, v[0]+1] )
        return self._element_constructor_( [v[0]+1, v[0]-1] )

    def _element_constructor_(self, *args, **keywords):
        return self.element_class(self,*args, **keywords)

    Element = IntegerPair




#################################################################
# Enumerated set of all kary trees
#################################################################
class KaryTrees_all(DisjointUnionEnumeratedSets, KaryTrees):

    def __init__(self):
        """
        TODO

        TESTS::

            sage: from sage.combinat.kary_tree import KaryTrees_all
            sage: K = KaryTrees_all()
            sage: K.cardinality()
            +Infinity

#            sage: it = iter(K)
#            sage: (next(it), next(it), next(it), next(it), next(it))
#            (., [., .], [., [., .]], [[., .], .], [., [., [., .]]])
#            sage: next(it).parent()
#            k-ary trees
            sage: K([None, None, None])
            [., ., .]

            TOTO : est-ce que l'on permet : K([]) et K(None) qui seraient égal
            respectivement à : K([None, None, None]) et K(None, 3).

            sage: K is KaryTrees_all()
            True
#            sage: TestSuite(K).run() # long time
            """

        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                IntegerPairs(), 
                lambda x : KaryTrees_size(arity=x[0], size=x[1])
            ),
            facade=True, keepkey = False
        )

    def _repr_(self):
        """
        TEST::

            sage: KaryTrees()   # indirect doctest
            k-ary trees
        """
        return "k-ary trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: K = KaryTrees()
            sage: 1 in K
            False
            sage: K([None, None, None]) in K
            True
        """
        return isinstance(x, self.element_class)

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TODO

        TESTS::

            sage: K = KaryTrees()
            sage: K([None, None])
            [., .]
        """
        return super(KaryTrees, self).__call__(x, *args, **keywords)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: KaryTrees().unlabelled_trees()
            3-ary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: KaryTrees().labelled_trees()
            Labelled k-ary trees
        """
        return LabelledKaryTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: K = KaryTrees(3)
            sage: K._element_constructor_([None, None, None])
            [., ., .]
            sage: K([[None, None, None],None, None]) # indirect doctest
            [[., ., .], ., .]
            sage: K(None)    # indirect doctest
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = KaryTree


#################################################################
# Enumerated set of all kary trees with a ficd arity
#################################################################
class KaryTrees_arity(DisjointUnionEnumeratedSets, KaryTrees):

    def __init__(self, arity):
        """
        TODO

        TESTS::

            sage: from sage.combinat.kary_tree import KaryTrees_all
            sage: K = KaryTrees_arity(3)
            sage: K.cardinality()
            +Infinity

#            sage: it = iter(K)
#            sage: (next(it), next(it), next(it), next(it), next(it))
#            (., [., .], [., [., .]], [[., .], .], [., [., [., .]]])
#            sage: next(it).parent()
#            k-ary trees
            sage: K([None, None, None])
            [., ., .]

            TOTO : est-ce que l'on permet : K([]) et K(None) qui seraient égal
            respectivement à : K([None, None, None]) et K(None, 3).

            sage: K is KaryTrees_arity(3)
            True
#            sage: TestSuite(K).run() # long time
            """
        self._arity = arity
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(), lambda x: KaryTrees_size(self._arity, x)
            ), facade=True, keepkey = False
        )

    def _repr_(self):
        """
        TEST::

            sage: KaryTrees(3)   # indirect doctest
            3-ary trees
        """
        return "%d-ary trees"%(self._arity)

    def __contains__(self, x):
        """
        TESTS::

            sage: K = KaryTrees(3)
            sage: 1 in K
            False
            sage: K([None, None, None]) in K
            True
            sage: K([None, None]) in K
            False
        """
        return isinstance(x, self.element_class) and ( self.is_empty() or self.arity() == self._arity )

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TODO

        TESTS::

            sage: K = KaryTrees()
            sage: K([None, None])
            [., .]
        """
        return super(KaryTrees, self).__call__(x, *args, **keywords)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: KaryTrees(3).unlabelled_trees()
            3-ary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: KaryTrees(3).labelled_trees()
            Labelled k-ary trees
        """
        return LabelledKaryTrees(self._arity)

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: K = KaryTrees(3)
            sage: K._element_constructor_([None, None, None])
            [., ., .]
            sage: K([[None, None, None],None, None]) # indirect doctest
            [[., ., .], ., .]
            sage: K(None)    # indirect doctest
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = KaryTree



from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
#################################################################
# Enumerated set of k-ary trees of a given size
#################################################################
class KaryTrees_size(KaryTrees):
    """
    The enumerated sets of k-ary trees of given size

    TESTS::

    TODO 

#        sage: from sage.combinat.kary_tree import KaryTrees_size
#        sage: for i in range(6): TestSuite(KaryTrees_size(3,i)).run()
    """
    def __init__(self, arity, size):
        """
        TESTS::

            sage: S = KaryTrees(2, 3)
            sage: S == loads(dumps(S))
            True

            sage: S is KaryTrees(2, 3)
            True
        """
        super(KaryTrees_size, self).__init__(
            category = FiniteEnumeratedSets()
        )
        self._arity = arity
        self._size = size

    def _repr_(self):
        """
        TESTS::

            sage: KaryTrees(2, 3)   # indirect doctest
            2-ary trees of size 3
        """
        return "%d-ary trees of size %s" % (self._arity, self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: K = KaryTrees(2, 3)
            sage: 1 in S
            False
            sage: K([[None, None],[None, None]]) in K
            True
        """
        return isinstance(x, self.element_class) and (
            x.is_empty() or (
                x.node_number() == self._size and x.arity() == self._arity
            )
        )

    def _an_element_(self):
        """
        TESTS::

            sage: KaryTrees(2, 0).an_element()  # indirect doctest
            .
        """
        return self.first()

    def cardinality(self):
        """
        The cardinality of ``self``

        TODO
        This is a ?? number.

        TESTS::

            TODO
            sage: KaryTrees(2, 0).cardinality()
            1
            sage: KaryTrees(2, 5).cardinality()
            42
        """
        return binomial(self._arity*self._size, self._size)/(
            (self._arity-1)*self._size + 1
        )

#    def random_element(self):
#        r"""
#        Return a random ``BinaryTree`` with uniform probability.
#
#        This method generates a random ``DyckWord`` and then uses a
#        bijection between Dyck words and binary trees.
#
#        EXAMPLES::
#
#            sage: BinaryTrees(5).random_element() # random
#            [., [., [., [., [., .]]]]]
#            sage: BinaryTrees(0).random_element()
#            .
#            sage: BinaryTrees(1).random_element()
#            [., .]
#
#        TESTS::
#
#            sage: all([BinaryTrees(10).random_element() in BinaryTrees(10) for i in range(20)])
#            True
#        """
#        from sage.combinat.dyck_word import CompleteDyckWords_size
#        return CompleteDyckWords_size(self._size).random_element().to_binary_tree()

    def __iter__(self):
        """
        A basic generator.

        .. TODO:: could be optimized.

        TESTS::

            sage: KaryTrees(2, 0).list()
            [.]
            sage: KaryTrees(2, 1).list()
            [[., .]]
            sage: KaryTrees(2, 3).list()
            [[., [., [., .]]], [., [[., .], .]], [[., .], [., .]], [[., [., .]], .], [[[., .], .], .]]
        """
        if self._size == 0:
            yield self._element_constructor_(None)
        else:
            for v in IntegerVectors( self._size-1, length=self._arity ):
                cp = CartesianProduct(
                    * map( lambda x:self.__class__(self._arity, x), v )
                )
                for children in cp:
                    yield self._element_constructor_( children )

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the elements generated by ``self``.

        TESTS::

            sage: S = KaryTrees(2, 3)
            sage: S._parent_for
            k-ary trees
        """
        return KaryTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: K = KaryTrees(2, 3)
            sage: K.element_class
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>
            sage: K.first().__class__ == KaryTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = KaryTrees(2, 0)
            sage: S([None])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
            sage: S(None)   # indirect doctest
            .

            sage: S = KaryTrees(2, 1)   # indirect doctest
            sage: S([None, None])
            [., .]
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.node_number() != self._size:
            raise ValueError("wrong number of nodes")
        if ( not res.is_empty() ) and ( res.arity() != self._arity ) :
            raise ValueError("wrong arity")
        return res

class LabelledKaryTree(AbstractLabelledClonableTree, KaryTree):
    pass
#    """
#    Labelled binary trees.
#
#    A labelled binary tree is a binary tree (see :class:`BinaryTree` for
#    the meaning of this) with a label assigned to each node.
#    The labels need not be integers, nor are they required to be distinct.
#    ``None`` can be used as a label.
#
#    .. WARNING::
#
#        While it is possible to assign values to leaves (not just nodes)
#        using this class, these labels are disregarded by various
#        methods such as
#        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.labels`,
#        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.map_labels`,
#        and (ironically)
#        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.leaf_labels`.
#
#    INPUT:
#
#    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
#      length `2` of labelled binary trees or convertible objects. This
#      corresponds to the standard recursive definition of a labelled
#      binary tree as being either a leaf, or a pair of:
#
#      - a pair of labelled binary trees,
#      - and a label.
#
#      (The label is specified in the keyword variable ``label``; see
#      below.)
#
#      Syntactic sugar allows leaving out all but the outermost calls
#      of the ``LabelledBinaryTree()`` constructor, so that, e. g.,
#      ``LabelledBinaryTree([LabelledBinaryTree(None),LabelledBinaryTree(None)])``
#      can be shortened to ``LabelledBinaryTree([None,None])``. However,
#      using this shorthand, it is impossible to label any vertex of
#      the tree other than the root (because there is no way to pass a
#      ``label`` variable without calling ``LabelledBinaryTree``
#      explicitly).
#
#      It is also allowed to abbreviate ``[None, None]`` by ``[]`` if
#      one does not want to label the leaves (which one should not do
#      anyway!).
#
#    - ``label`` -- (default: ``None``) the label to be put on the root
#      of this tree.
#
#    - ``check`` -- (default: ``True``) whether checks should be
#      performed or not.
#
#    .. TODO::
#
#        It is currently not possible to use ``LabelledBinaryTree()``
#        as a shorthand for ``LabelledBinaryTree(None)`` (in analogy to
#        similar syntax in the ``BinaryTree`` class).
#
#    EXAMPLES::
#
#        sage: LabelledBinaryTree(None)
#        .
#        sage: LabelledBinaryTree(None, label="ae")    # not well supported
#        'ae'
#        sage: LabelledBinaryTree([])
#        None[., .]
#        sage: LabelledBinaryTree([], label=3)    # not well supported
#        3[., .]
#        sage: LabelledBinaryTree([None, None])
#        None[., .]
#        sage: LabelledBinaryTree([None, None], label=5)
#        5[., .]
#        sage: LabelledBinaryTree([None, []])
#        None[., None[., .]]
#        sage: LabelledBinaryTree([None, []], label=4)
#        4[., None[., .]]
#        sage: LabelledBinaryTree([[], None])
#        None[None[., .], .]
#        sage: LabelledBinaryTree("[[], .]", label=False)
#        False[None[., .], .]
#        sage: LabelledBinaryTree([None, LabelledBinaryTree([None, None], label=4)], label=3)
#        3[., 4[., .]]
#        sage: LabelledBinaryTree([None, BinaryTree([None, None])], label=3)
#        3[., None[., .]]
#
#        sage: LabelledBinaryTree([[], None, []])
#        Traceback (most recent call last):
#        ...
#        ValueError: this is not a binary tree
#
#        sage: LBT = LabelledBinaryTree
#        sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
#        4[None[2[., .], .], .]
#
#    TESTS::
#
#        sage: t1 = LabelledBinaryTree([[None, [[],[[], None]]],[[],[]]])
#        sage: t2 = LabelledBinaryTree([[[],[]],[]])
#        sage: with t1.clone() as t1c:
#        ....:     t1c[1,1,1] = t2
#        sage: t1 == t1c
#        False
#
#    We check for :trac:`16314`::
#
#        sage: t1 = LBT([ LBT([LBT([], label=2),
#        ....:                 LBT([], label=5)], label=6),
#        ....:            None], label=4); t1
#        4[6[2[., .], 5[., .]], .]
#        sage: class Foo(LabelledBinaryTree):
#        ....:     pass
#        sage: t2 = Foo(t1.parent(), t1); t2
#        4[6[2[., .], 5[., .]], .]
#        sage: t2.label()
#        4
#        sage: t2[0].label()
#        6
#        sage: t2.__class__, t2[0].__class__
#        (<class '__main__.Foo'>, <class '__main__.Foo'>)
#    """
#    @staticmethod
#    def __classcall_private__(cls, *args, **opts):
#        """
#        Ensure that trees created by the sets and directly are the same and
#        that they are instances of :class:`LabelledTree`.
#
#        TESTS::
#
#            sage: issubclass(LabelledBinaryTrees().element_class, LabelledBinaryTree)
#            True
#            sage: t0 = LabelledBinaryTree([[],[[], None]], label = 3)
#            sage: t0.parent()
#            Labelled binary trees
#            sage: type(t0)
#            <class 'sage.combinat.binary_tree.LabelledBinaryTrees_with_category.element_class'>
#        """
#        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)
#
#    @lazy_class_attribute
#    def _auto_parent(cls):
#        """
#        The automatic parent of the elements of this class.
#
#        When calling the constructor of an element of this class, one needs a
#        parent. This class attribute specifies which parent is used.
#
#        EXAMPLES::
#
#            sage: LabelledBinaryTree._auto_parent
#            Labelled binary trees
#            sage: LabelledBinaryTree([], label = 3).parent()
#            Labelled binary trees
#        """
#        return LabelledBinaryTrees()
#
#    def _repr_(self):
#        """
#        TESTS::
#
#            sage: LBT = LabelledBinaryTree
#            sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
#            4[None[2[., .], .], .]
#            sage: LBT([[],[[], None]], label = 3)   # indirect doctest
#            3[None[., .], None[None[., .], .]]
#        """
#        if not self:
#            if self._label is not None:
#                return repr(self._label)
#            else:
#                return "."
#        else:
#            return "%s%s"%(self._label, self[:])
#
#    def binary_search_insert(self, letter):
#        """
#        Return the result of inserting a letter ``letter`` into the
#        right strict binary search tree ``self``.
#
#        INPUT:
#
#        - ``letter`` -- any object comparable with the labels of ``self``
#
#        OUTPUT:
#
#        The right strict binary search tree ``self`` with ``letter``
#        inserted into it according to the binary search insertion
#        algorithm.
#
#        .. NOTE:: ``self`` is supposed to be a binary search tree.
#                  This is not being checked!
#
#        A right strict binary search tree is defined to be a labelled
#        binary tree such that for each node `n` with label `x`,
#        every descendant of the left child of `n` has a label `\leq x`,
#        and every descendant of the right child of `n` has a label
#        `> x`. (Here, only nodes count as descendants, and every node
#        counts as its own descendant too.) Leaves are assumed to have
#        no labels.
#
#        Given a right strict binary search tree `t` and a letter `i`,
#        the result of inserting `i` into `t` (denoted `Ins(i, t)` in
#        the following) is defined recursively as follows:
#
#        - If `t` is empty, then `Ins(i, t)` is the tree with one node
#          only, and this node is labelled with `i`.
#
#        - Otherwise, let `j` be the label of the root of `t`. If
#          `i > j`, then `Ins(i, t)` is obtained by replacing the
#          right child of `t` by `Ins(i, r)` in `t`, where `r` denotes
#          the right child of `t`. If `i \leq j`, then `Ins(i, t)` is
#          obtained by replacing the left child of `t` by `Ins(i, l)`
#          in `t`, where `l` denotes the left child of `t`.
#
#        See, for example, [HNT05]_ for properties of this algorithm.
#
#        .. WARNING::
#
#            If `t` is nonempty, then inserting `i` into `t` does not
#            change the root label of `t`. Hence, as opposed to
#            algorithms like Robinson-Schensted-Knuth, binary
#            search tree insertion involves no bumping.
#
#        EXAMPLES:
#
#        The example from Fig. 2 of [HNT05]_::
#
#            sage: LBT = LabelledBinaryTree
#            sage: x = LBT(None)
#            sage: x
#            .
#            sage: x = x.binary_search_insert("b"); x
#            b[., .]
#            sage: x = x.binary_search_insert("d"); x
#            b[., d[., .]]
#            sage: x = x.binary_search_insert("e"); x
#            b[., d[., e[., .]]]
#            sage: x = x.binary_search_insert("a"); x
#            b[a[., .], d[., e[., .]]]
#            sage: x = x.binary_search_insert("b"); x
#            b[a[., b[., .]], d[., e[., .]]]
#            sage: x = x.binary_search_insert("d"); x
#            b[a[., b[., .]], d[d[., .], e[., .]]]
#            sage: x = x.binary_search_insert("a"); x
#            b[a[a[., .], b[., .]], d[d[., .], e[., .]]]
#            sage: x = x.binary_search_insert("c"); x
#            b[a[a[., .], b[., .]], d[d[c[., .], .], e[., .]]]
#
#        Other examples::
#
#            sage: LBT = LabelledBinaryTree
#            sage: LBT(None).binary_search_insert(3)
#            3[., .]
#            sage: LBT([], label = 1).binary_search_insert(3)
#            1[., 3[., .]]
#            sage: LBT([], label = 3).binary_search_insert(1)
#            3[1[., .], .]
#            sage: res = LBT(None)
#            sage: for i in [3,1,5,2,4,6]:
#            ....:     res = res.binary_search_insert(i)
#            sage: res
#            3[1[., 2[., .]], 5[4[., .], 6[., .]]]
#        """
#        LT = self.parent()._element_constructor_
#        if not self:
#            return LT([], label = letter)
#        else:
#            if letter <= self.label():
#                fils = self[0].binary_search_insert(letter)
#                return LT([fils, self[1]], label=self.label())
#            else:
#                fils = self[1].binary_search_insert(letter)
#                return LT([self[0], fils], label=self.label())
#
#    def semistandard_insert(self, letter):
#        """
#        Return the result of inserting a letter ``letter`` into the
#        semistandard tree ``self`` using the bumping algorithm.
#
#        INPUT:
#
#        - ``letter`` -- any object comparable with the labels of ``self``
#
#        OUTPUT:
#
#        The semistandard tree ``self`` with ``letter`` inserted into it
#        according to the bumping algorithm.
#
#        .. NOTE:: ``self`` is supposed to be a semistandard tree.
#                  This is not being checked!
#
#        A semistandard tree is defined to be a labelled binary tree
#        such that for each node `n` with label `x`, every descendant of
#        the left child of `n` has a label `> x`, and every descendant
#        of the right child of `n` has a label `\geq x`. (Here, only
#        nodes count as descendants, and every node counts as its own
#        descendant too.) Leaves are assumed to have no labels.
#
#        Given a semistandard tree `t` and a letter `i`, the result of
#        inserting `i` into `t` (denoted `Ins(i, t)` in the following)
#        is defined recursively as follows:
#
#        - If `t` is empty, then `Ins(i, t)` is the tree with one node
#          only, and this node is labelled with `i`.
#
#        - Otherwise, let `j` be the label of the root of `t`. If
#          `i \geq j`, then `Ins(i, t)` is obtained by replacing the
#          right child of `t` by `Ins(i, r)` in `t`, where `r` denotes
#          the right child of `t`. If `i < j`, then `Ins(i, t)` is
#          obtained by replacing the label at the root of `t` by `i`,
#          and replacing the left child of `t` by `Ins(j, l)`
#          in `t`, where `l` denotes the left child of `t`.
#
#        This algorithm is similar to the Robinson-Schensted-Knuth
#        insertion algorithm for semistandard Young tableaux.
#
#        AUTHORS:
#
#        - Darij Grinberg (10 Nov 2013).
#
#        EXAMPLES::
#
#            sage: LBT = LabelledBinaryTree
#            sage: x = LBT(None)
#            sage: x
#            .
#            sage: x = x.semistandard_insert("b"); x
#            b[., .]
#            sage: x = x.semistandard_insert("d"); x
#            b[., d[., .]]
#            sage: x = x.semistandard_insert("e"); x
#            b[., d[., e[., .]]]
#            sage: x = x.semistandard_insert("a"); x
#            a[b[., .], d[., e[., .]]]
#            sage: x = x.semistandard_insert("b"); x
#            a[b[., .], b[d[., .], e[., .]]]
#            sage: x = x.semistandard_insert("d"); x
#            a[b[., .], b[d[., .], d[e[., .], .]]]
#            sage: x = x.semistandard_insert("a"); x
#            a[b[., .], a[b[d[., .], .], d[e[., .], .]]]
#            sage: x = x.semistandard_insert("c"); x
#            a[b[., .], a[b[d[., .], .], c[d[e[., .], .], .]]]
#
#        Other examples::
#
#            sage: LBT = LabelledBinaryTree
#            sage: LBT(None).semistandard_insert(3)
#            3[., .]
#            sage: LBT([], label = 1).semistandard_insert(3)
#            1[., 3[., .]]
#            sage: LBT([], label = 3).semistandard_insert(1)
#            1[3[., .], .]
#            sage: res = LBT(None)
#            sage: for i in [3,1,5,2,4,6]:
#            ....:     res = res.semistandard_insert(i)
#            sage: res
#            1[3[., .], 2[5[., .], 4[., 6[., .]]]]
#        """
#        LT = self.parent()._element_constructor_
#        if not self:
#            return LT([], label = letter)
#        else:
#            root_label = self.label()
#            if letter < root_label:
#                fils = self[0].semistandard_insert(root_label)
#                return LT([fils, self[1]], label=letter)
#            else:
#                fils = self[1].semistandard_insert(letter)
#                return LT([self[0], fils], label=root_label)
#
#    def right_rotate(self):
#        r"""
#        Return the result of right rotation applied to the labelled
#        binary tree ``self``.
#
#        Right rotation on labelled binary trees is defined as
#        follows: Let `T` be a labelled binary tree such that the
#        left child of the root of `T` is a node. Let
#        `C` be the right child of the root of `T`, and let `A`
#        and `B` be the left and right children of the left child
#        of the root of `T`. (Keep in mind that nodes of trees are
#        identified with the subtrees consisting of their
#        descendants.) Furthermore, let `y` be the label at the
#        root of `T`, and `x` be the label at the left child of the
#        root of `T`.
#        Then, the right rotation of `T` is the labelled binary
#        tree in which the root is labelled `x`, the left child of
#        the root is `A`, whereas the right child of the root is a
#        node labelled `y` whose left and right children are `B`
#        and `C`. In pictures::
#
#            |     y                      x     |
#            |    / \                    / \    |
#            |   x   C -right-rotate->  A   y   |
#            |  / \                        / \  |
#            | A   B                      B   C |
#
#        Right rotation is the inverse operation to left rotation
#        (:meth:`left_rotate`).
#
#        TESTS::
#
#            sage: LB = LabelledBinaryTree
#            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
#            y[x[A[., .], B[., .]], C[., .]]
#            sage: b.right_rotate()
#            x[A[., .], y[B[., .], C[., .]]]
#        """
#        B = self.parent()._element_constructor_
#        s0 = self[0]
#        return B([s0[0], B([s0[1], self[1]], self.label())], s0.label())
#
#    def left_rotate(self):
#        r"""
#        Return the result of left rotation applied to the labelled
#        binary tree ``self``.
#
#        Left rotation on labelled binary trees is defined as
#        follows: Let `T` be a labelled binary tree such that the
#        right child of the root of `T` is a node. Let
#        `A` be the left child of the root of `T`, and let `B`
#        and `C` be the left and right children of the right child
#        of the root of `T`. (Keep in mind that nodes of trees are
#        identified with the subtrees consisting of their
#        descendants.) Furthermore, let `x` be the label at the
#        root of `T`, and `y` be the label at the right child of the
#        root of `T`.
#        Then, the left rotation of `T` is the labelled binary tree
#        in which the root is labelled `y`, the right child of the
#        root is `C`, whereas the left child of the root is a node
#        labelled `x` whose left and right children are `A` and `B`.
#        In pictures::
#
#           |     y                    x     |
#           |    / \                  / \    |
#           |   x   C <-left-rotate- A   y   |
#           |  / \                      / \  |
#           | A   B                    B   C |
#
#        Left rotation is the inverse operation to right rotation
#        (:meth:`right_rotate`).
#
#        TESTS::
#
#            sage: LB = LabelledBinaryTree
#            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
#            y[x[A[., .], B[., .]], C[., .]]
#            sage: b == b.right_rotate().left_rotate()
#            True
#        """
#        B = self.parent()._element_constructor_
#        s1 = self[1]
#        return B([B([self[0], s1[0]], self.label()), s1[1]], s1.label())
#
#    def heap_insert(self, l):
#        r"""
#        Return the result of inserting a letter ``l`` into the binary
#        heap (tree) ``self``.
#
#        A binary heap is a labelled complete binary tree such that for
#        each node, the label at the node is greater or equal to the
#        label of each of its child nodes. (More precisely, this is
#        called a max-heap.)
#
#        For example::
#
#            |     _7_   |
#            |    /   \  |
#            |   5     6 |
#            |  / \      |
#            | 3   4     |
#
#        is a binary heap.
#
#        See :wikipedia:`Binary_heap#Insert` for a description of how to
#        insert a letter into a binary heap. The result is another binary
#        heap.
#
#        INPUT:
#
#        - ``letter`` -- any object comparable with the labels of ``self``
#
#        .. NOTE::
#
#            ``self`` is assumed to be a binary heap (tree). No check is
#            performed.
#
#        TESTS::
#
#            sage: h = LabelledBinaryTree(None)
#            sage: h = h.heap_insert(3); ascii_art([h])
#            [ 3 ]
#            sage: h = h.heap_insert(4); ascii_art([h])
#            [   4 ]
#            [  /  ]
#            [ 3   ]
#            sage: h = h.heap_insert(6); ascii_art([h])
#            [   6   ]
#            [  / \  ]
#            [ 3   4 ]
#            sage: h = h.heap_insert(2); ascii_art([h])
#            [     6   ]
#            [    / \  ]
#            [   3   4 ]
#            [  /      ]
#            [ 2       ]
#            sage: ascii_art([h.heap_insert(5)])
#            [     _6_   ]
#            [    /   \  ]
#            [   5     4 ]
#            [  / \      ]
#            [ 2   3     ]
#        """
#        B = self.parent()._element_constructor_
#        if self.is_empty():
#            return B([], l)
#
#        if self.label() < l:
#            label_root = l
#            label_insert = self.label()
#        else:
#            label_root = self.label()
#            label_insert = l
#        L, R = self
#        dL = L.depth()
#        dR = R.depth()
#        # if depth of L is greater than the depth of R
#        if dL > dR:
#            # if L is perfect we insert in R
#            if L.is_perfect():
#                return B([L, R.heap_insert(label_insert)], label_root)
#            # we insert in L
#            return B([L.heap_insert(label_insert), R], label_root)
#        # else ==> dL == dR
#        # if R is perfect we have to insert on the leftmost leaf
#        if R.is_perfect():
#            # ## TODO:: can be optimized...
#            return B([L.heap_insert(label_insert), R], label_root)
#        # else we insert on the right
#        return B([L, R.heap_insert(label_insert)], label_root)
#
#    _UnLabelled = BinaryTree


class LabelledKaryTrees(LabelledOrderedTrees):
    pass
#    """
#    This is a parent stub to serve as a factory class for trees with various
#    labels constraints.
#    """
#    def _repr_(self):
#        """
#        TESTS::
#
#            sage: LabelledBinaryTrees()   # indirect doctest
#            Labelled binary trees
#        """
#        return "Labelled binary trees"
#
#    def _an_element_(self):
#        """
#        Return a labelled binary tree.
#
#        EXAMPLE::
#
#            sage: LabelledBinaryTrees().an_element()   # indirect doctest
#            toto[42[3[., .], 3[., .]], 5[None[., .], None[., .]]]
#        """
#        LT = self._element_constructor_
#        t  = LT([], label = 3)
#        t1 = LT([t,t], label = 42)
#        t2  = LT([[], []], label = 5)
#        return LT([t1,t2], label = "toto")
#
#    def unlabelled_trees(self):
#        """
#        Return the set of unlabelled trees associated to ``self``.
#
#        EXAMPLES::
#
#            sage: LabelledBinaryTrees().unlabelled_trees()
#            Binary trees
#
#        This is used to compute the shape::
#
#            sage: t = LabelledBinaryTrees().an_element().shape(); t
#            [[[., .], [., .]], [[., .], [., .]]]
#            sage: t.parent()
#            Binary trees
#
#        TESTS::
#
#            sage: t = LabelledBinaryTrees().an_element()
#            sage: t.canonical_labelling()
#            4[2[1[., .], 3[., .]], 6[5[., .], 7[., .]]]
#        """
#        return BinaryTrees_all()
#
#    def labelled_trees(self):
#        """
#        Return the set of labelled trees associated to ``self``.
#
#        EXAMPLES::
#
#            sage: LabelledBinaryTrees().labelled_trees()
#            Labelled binary trees
#        """
#        return self
#
#    Element = LabelledBinaryTree




