r"""
TESTS::

    sage: from sage.combinat.words.word import FiniteWord_char
    sage: W = Words(IntegerRange(0,256))
    sage: w = FiniteWord_char(W, range(1,10)*2)
    sage: w
    word: 123456789123456789

    sage: w == w
    True
    sage: w != w
    False
    sage: w[:-1] != w[1:]
    True
    sage: w < w[1:] and w[1:] > w
    True
    sage: w > w[1:] or w[1:] < w
    False

    sage: list(w) == [w[i] for i in range(len(w))]
    True

    sage: type(len(w))
    <type 'int'>
    sage: type(w.length())
    <type 'sage.rings.integer.Integer'>

    sage: w.has_prefix([1,2,3,4])
    True
    sage: w.has_prefix([1,2,4,4])
    False
    sage: w.has_prefix(FiniteWord_char(W,[1,2,3,4]))
    True
    sage: w.has_prefix(FiniteWord_char(W,[1,2,4,4]))
    False

    sage: w.is_palindrome()
    False
    sage: (w*w[::-1]).is_palindrome()
    True
    sage: (w[:-1:]*w[::-1]).is_palindrome()
    True

    sage: w.is_lyndon()
    False
    sage: FiniteWord_char(W, range(10)+[10,10]).is_lyndon()
    True

    sage: w.is_square()
    True
    sage: w[:-2].is_square()
    False
    sage: w.is_square_free()
    False
    sage: w[:-1].is_square_free()
    True
    sage: u = FiniteWord_char(W, [randint(0,255) for i in range(10)])
    sage: (u*u).is_square()
    True
    sage: (u*u*u).is_cube()
    True

    sage: len(w.factor_set())
    127
    sage: w.rauzy_graph(5)
    Looped digraph on 9 vertices

    sage: u = FiniteWord_char(W,[1,2,3])
    sage: u.first_pos_in(w)
    0
    sage: u.first_pos_in(w[1:])
    8

    sage: w = FiniteWord_char(W, [0,1,2,3])
    sage: w
    word: 0123
    sage: w ** (1/2)
    word: 01
    sage: w ** 2
    word: 01230123
    sage: w ** 3
    word: 012301230123
    sage: w ** (7/2)
    word: 01230123012301
    sage: len(((w ** 2) ** 3) ** 5) == len(w) * 2 * 3 * 5
    True
"""

include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

cimport cython
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.rational cimport Rational
from libc.string cimport memcpy, memcmp

cdef extern from "Python.h":
    # check functions
    int PyIndex_Check(object o)
    int PySlice_Check(object o)
    int PySequence_Check(object o)
    int PyNumber_Check(object o)

    # get numbers from Python slice
    int PySlice_GetIndicesEx(object slice, Py_ssize_t length,
            Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
            Py_ssize_t *slicelength)

# the maximum value of a size_t
cdef size_t SIZE_T_MAX = -(<size_t> 1)

def reversed_word_iterator(Word_char w):
    cdef ssize_t i
    for i in range(w._length-1, 0, -1):
        yield w._data[i]
    yield w._data[0]

cdef class Word_char(object):
    r"""
    A Fast class for words.

    Currently, only handles letters in [0,256].
    """
    cdef unsigned char * _data
    cdef size_t _length
    cdef Word_char _master
    cdef public object _parent
    cdef object _hash

    def __cinit__(self):
        self._master = None
        self._hash = None
        self._data = NULL
        self._length = 0

    def __init__(self, parent, data):
        self._parent = parent

        if data:
            if not PySequence_Check(data):
                raise TypeError("not able to initialize a word from {}".format(data))

            self._set_data(data)

    @cython.boundscheck(False) # assume that indexing will not cause any IndexErrors
    @cython.wraparound(False)  # not check not correctly handle negative indices
    cdef _set_data(self, data):
        r"""
        set the attribute ._data and ._length from the sequence data
        (usually data is a word, a tuple or a list)
        """
        cdef size_t i
        self._length = len(data)
        self._data = <unsigned char *> sage_malloc(self._length * sizeof(unsigned char))
        if self._data == NULL:
            raise MemoryError

        for i in range(self._length):
            self._data[i] = data[i]

    def __dealloc__(self):
        if self._master is None:
            sage_free(self._data)

    def __nonzero__(self):
        return self._length != 0

    def __len__(self):
        return self._length

    def length(self):
        return smallInteger(self._length)

    # TO DISCUSS: in Integer (sage.rings.integer) this method is actually an
    # external function. But we might want to have several possible inheritance.
    cdef _new_c(self, unsigned char * data, size_t length, Word_char master):
        cdef Word_char other = PY_NEW_SAME_TYPE(self)
        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__
        other._data = data
        other._master = master # can be None
        other._length = length
        other._parent = self._parent

        return other

    def __hash__(self):
        r"""
        Return the hash value.
        """
        cdef int res = 5381
        cdef size_t i
        if self._hash is None:
            for i in range(self._length):
                res = ((res << 5) + res) + self._data[i]
            self._hash = res
        return self._hash

    def __richcmp__(self, other, op):
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=
        if not PY_TYPE_CHECK(other, Word_char):
            return NotImplemented

        # word of different lengths are not equal!
        if (op == 2 or op == 3) and (<Word_char> self)._length != (<Word_char> other)._length:
            return op == 3

        cdef int test = (<Word_char> self)._lexico_cmp(other)
        if test < 0:
            return op < 2 or op == 3
        elif test > 0:
            return op > 3
        else:
            return op == 1 or op == 2 or op == 5

    def __cmp__(self, other):
        if not PY_TYPE_CHECK(other, Word_char):
            return NotImplemented

        cdef int test = self._lexico_cmp(other)
        if test:
            return test
        return (<Py_ssize_t> self._length) - (<Py_ssize_t> (<Word_char> other)._length)

    cdef int _lexico_cmp(self, Word_char other) except -2:
        r"""
        Lexicographic comparison of self and other up to
        the letter at position min(len(self),len(other))
        """
        cdef size_t l = min(self._length, other._length)

        sig_on()
        cdef int test = memcmp(<void *> (<Word_char> self)._data,
                      <void *> (<Word_char> other)._data,
                      l * sizeof(unsigned char))
        sig_off()

        if test == 0:
            return 0
        if test < 0:
            return -1
        else:
            return 1

    def _get_info(self):
        r"""
        Temporary function that print useful debug information
        """
        if self._master:
            print "master at %u"%(id(self._master))
        print "data at %u"%(<size_t>self._data)

    def __getitem__(self, key):
        cdef Py_ssize_t i, start, stop, step, slicelength
        cdef unsigned char * data
        cdef size_t j,k
        if PySlice_Check(key):
            # here the key is a slice
            if PySlice_GetIndicesEx(key,
                    self._length,
                    &start, &stop, &step,
                    &slicelength) < 0:
                return  # this automatically raise an error because
                        # PySlice_GetIndices already did the job
            if step == 1:
                return self._new_c(self._data+start, stop-start, self)
            data = <unsigned char *> sage_malloc(slicelength * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[j] = self._data[k]
                j += 1
            return self._new_c(data, slicelength, None)

        elif PyIndex_Check(key):
            # here the key is an int
            i = key
            if i < 0:
                i += self._length;
            if i < 0 or i >= self._length:
                raise IndexError("word index out of range")
            return self._data[i]

        raise TypeError("word indices must be integers")

    def __iter__(self):
        cdef size_t i
        for i in range(self._length):
            yield self._data[i]

    def __reversed__(self):
        return reversed_word_iterator(self)

    cdef _concatenate(self, Word_char other):
        cdef unsigned char * data
        data = <unsigned char *> sage_malloc((self._length + other._length) * sizeof(unsigned char))
        if data == NULL:
            raise MemoryError

        memcpy(data, self._data, self._length * sizeof(unsigned char))
        memcpy(data+self._length, other._data, other._length * sizeof(unsigned char))

        return self._new_c(data, self._length + other._length, None)

    def __mul__(self, other):
        r"""
        Concatenation of ``self`` and ``other``.

        The result is automatically converted to a Word_char.
        """
        cdef Word_char w

        if PY_TYPE_CHECK(other, Word_char):
            return (<Word_char> self)._concatenate(other)

        elif PySequence_Check(other):
            # we convert other to a Word_char and perform the concatenation
            w = (<Word_char> self)._new_c(NULL, 0, None)
            w._set_data(other)
            return (<Word_char> self)._concatenate(w)

        raise TypeError("not able to initialize a word from {}".format(other))

    def __pow__(self, exp, mod):
        r"""
        Power
        """
        if not PyNumber_Check(exp):
            raise ValueError("the exponent must be a number or infinity")
        if mod is not None:
            raise ValueError("a word can not be taken modulo")

        if exp == float('inf'):
            from sage.rings.infinity import Infinity
            fcn = lambda n: self[n % self.length()]
            return self._parent(fcn, length=Infinity)

        if exp < 0:
            raise ValueError("can not take negative power of a word")

        cdef Word_char w = self
        cdef size_t i, rest

        if PY_TYPE_CHECK_EXACT(exp, Rational):
            if w._length % exp.denominator():
                raise ValueError("undefined")
            i = exp.floor()
            rest = (exp - exp.floor()) * w._length
        else:
            i = exp
            rest = 0

        # first handle the cases (i*length + rest) <= length and return the
        # corresponding prefix of self
        if i == 1 and rest == 0:
            return self
        if w._length == 0:
            return w._new_c(NULL, 0, None)
        if i == 0:
            if rest == 0:
                return w._new_c(NULL, 0, None)
            else:
                return w._new_c(w._data, rest, self)

        # now consider non trivial powers
        if w._length > SIZE_T_MAX / (i+1):
            raise OverflowError("the length of the result is too large")
        cdef size_t new_length = w._length * i + rest
        cdef unsigned char * data = <unsigned char *> sage_malloc(new_length * sizeof(unsigned char))
        if data == NULL:
            raise MemoryError

        cdef Py_ssize_t j = w._length
        memcpy(data, w._data, j * sizeof(unsigned char))
        while 2*j < new_length:
            memcpy(data + j, data, j * sizeof(unsigned char))
            j *= 2
        memcpy(data + j, data, (new_length - j) * sizeof(unsigned char))

        return w._new_c(data, new_length, None)

    @cython.boundscheck(False)
    def has_prefix(self, other):
        r"""
        Test whether ``other`` is a prefix of ``self``.

        INPUT:

        - ``other`` -- a word or a sequence (e.g. tuple, list)
        """
        cdef size_t i
        cdef int test
        cdef unsigned char * data
        cdef Word_char w

        if PY_TYPE_CHECK(other, Word_char):
            # C level
            w = <Word_char> other
            if w._length > self._length:
                return False
            return memcmp(self._data, w._data, w._length) == 0

        elif PySequence_Check(other):
            # python level
            if len(other) > self._length:
                return False

            for i in range(len(other)):
                if other[i] != self._data[i]:
                    return False
            return True

        raise TypeError("not able to initialize a word from {}".format(other))
