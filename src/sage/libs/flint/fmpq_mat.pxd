# distutils: libraries = flint
# distutils: depends = flint/fmpq_mat.h

from sage.libs.flint.types cimport fmpz_t, fmpz, fmpq_t, fmpq, fmpz_mat_t, fmpq_mat_t, flint_rand_t, mp_bitcnt_t

# cdef extern from "flint/fmpq_mat.h"
cdef extern from "flint_wrap.h":
    fmpq * fmpq_mat_entry(const fmpq_mat_t mat, long i, long j)
    fmpz * fmpq_mat_entry_num(const fmpq_mat_t mat, long i, long j)
    fmpz * fmpq_mat_entry_den(const fmpq_mat_t mat, long i, long j)
    long fmpq_mat_nrows(const fmpq_mat_t mat)
    long fmpq_mat_ncols(const fmpq_mat_t mat)
    void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
    void fmpq_mat_clear(fmpq_mat_t mat)
    void fmpq_mat_swap(fmpq_mat_t mat1, fmpq_mat_t mat2)
    void fmpq_mat_window_init(fmpq_mat_t window, const fmpq_mat_t mat, long r1, long c1, long r2, long c2)
    void fmpq_mat_window_clear(fmpq_mat_t window)
    void fmpq_mat_concat_horizontal(fmpq_mat_t res, const fmpq_mat_t mat1,  const fmpq_mat_t mat2)
    void fmpq_mat_concat_vertical(fmpq_mat_t res, const fmpq_mat_t mat1,  const fmpq_mat_t mat2)
    void fmpq_mat_print(const fmpq_mat_t mat)
    void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_hilbert_matrix(fmpq_mat_t mat)
    void fmpq_mat_set(fmpq_mat_t dest, const fmpq_mat_t src)
    void fmpq_mat_zero(fmpq_mat_t mat)
    void fmpq_mat_one(fmpq_mat_t mat)
    void fmpq_mat_transpose(fmpq_mat_t rop, const fmpq_mat_t op)
    void fmpq_mat_add(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
    void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
    void fmpq_mat_neg(fmpq_mat_t rop, const fmpq_mat_t op)
    void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)
    void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)
    bint fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2)
    bint fmpq_mat_is_integral(const fmpq_mat_t mat)
    bint fmpq_mat_is_zero(const fmpq_mat_t mat)
    bint fmpq_mat_is_empty(const fmpq_mat_t mat)
    bint fmpq_mat_is_square(const fmpq_mat_t mat)
    int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, const fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den, const fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, const fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2)
    void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, const fmpq_mat_t mat, const fmpz_t mod)
    void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src)
    void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t div)
    int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t mod)
    void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
    void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
    void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
    void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A, const fmpz_mat_t B)
    void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A, const fmpq_mat_t B)
    void fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat)
    void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat)
    int fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
    int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
    int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
    int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A)
    int fmpq_mat_pivot(long * perm, fmpq_mat_t mat, long r, long c)
    long fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A)
    long fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A)
    long fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)
    void fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A)

