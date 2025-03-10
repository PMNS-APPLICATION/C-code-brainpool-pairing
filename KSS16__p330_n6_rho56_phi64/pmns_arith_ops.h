#ifndef POLY_MULT_ADD
#define POLY_MULT_ADD


void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void add_poly(int64_t *rop, int64_t *pa, int64_t *pb);

void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void neg_poly(int64_t *rop, int64_t *op);
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar);
void double_poly_coeffs(int64_t *rop, int64_t *op);
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos);

void add_lpoly(int128_t *rop, int128_t *pa, int128_t *pb);

void scalar_mult_lpoly(int128_t *rop, int64_t *op, uint64_t scalar);

void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb);

void square_mod_poly(int64_t *rop, int64_t *pa);

void exact_coeffs_reduction(int64_t *rop, int64_t *op);

void from_mont_domain(int64_t *rop, int64_t *op);

int64_t equality_check_polys(int64_t *pA, int64_t *pB);

void internal_reduction(int64_t *rop, int128_t *op);

#endif

