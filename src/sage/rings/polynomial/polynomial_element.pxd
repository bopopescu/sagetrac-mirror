from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer
from .polynomial_compiled cimport CompiledPolynomialFunction


cdef class Polynomial(CommutativeAlgebraElement):
    cdef Polynomial _new_generic(self, list coeffs)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cpdef Polynomial truncate(self, long n)
    cpdef Polynomial inverse_series_trunc(self, long prec)
    cdef long _hash_c(self) except -1
    cpdef constant_coefficient(self)
    cpdef Polynomial _new_constant_poly(self, a, Parent P)
    cpdef list list(self, bint copy=*)
    cpdef _mul_generic(self, right)
    cdef _square_generic(self)

    cpdef bint is_zero(self)
    cpdef bint is_one(self)

    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, right)
    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n)
    cpdef Polynomial _power_trunc(self, unsigned long n, long prec)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

    cdef get_unsafe(self, Py_ssize_t i)
    cpdef long number_of_terms(self)

    cdef public dict __cached_methods

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P)
    cdef list __coeffs
    cdef int __normalize(self) except -1

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    pass

cpdef is_Polynomial(f)
cpdef Polynomial generic_power_trunc(Polynomial p, Integer n, long prec)
cpdef list _dict_to_list(dict x, zero)

