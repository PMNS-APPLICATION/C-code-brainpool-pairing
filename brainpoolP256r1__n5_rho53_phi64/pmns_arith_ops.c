#include "pmns_arith_ops.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

//~ --------------------------------------------------------------------

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ --------------------------------------------------------------------

//~ computes : pa + 2.pb
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + 2*pb[j];
}

//~ computes : pa - 2.pb
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - 2*pb[j];
}

//~ --------------------------------------------------------------------

void add_lpoly(int128_t *rop, int128_t *pa, int128_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_lpoly(int128_t *rop, int64_t *op, uint64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = (int128_t)op[j] * scalar;
}

//~ --------------------------------------------------------------------

//~ Computes: pa*pb mod(E)
//~ Note: E(X) = X^5 - 5
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pb[0] + (((int128_t)pa[1] * pb[4] + (int128_t)pa[2] * pb[3] + (int128_t)pa[3] * pb[2] + (int128_t)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128_t)pa[0] * pb[1] + (int128_t)pa[1] * pb[0] + (((int128_t)pa[2] * pb[4] + (int128_t)pa[3] * pb[3] + (int128_t)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128_t)pa[0] * pb[2] + (int128_t)pa[1] * pb[1] + (int128_t)pa[2] * pb[0] + (((int128_t)pa[3] * pb[4] + (int128_t)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128_t)pa[0] * pb[3] + (int128_t)pa[1] * pb[2] + (int128_t)pa[2] * pb[1] + (int128_t)pa[3] * pb[0] + (((int128_t)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128_t)pa[0] * pb[4] + (int128_t)pa[1] * pb[3] + (int128_t)pa[2] * pb[2] + (int128_t)pa[3] * pb[1] + (int128_t)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes: pa^2 mod(E)
//~ Note: E(X) = X^5 - 5
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pa[0] + (((int128_t)pa[3] * pa[2] + (int128_t)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128_t)pa[1] * pa[0]) << 1) + (((((int128_t)pa[4] * pa[2]) << 1) + (int128_t)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128_t)pa[2] * pa[0]) << 1) + (int128_t)pa[1] * pa[1] + (((int128_t)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128_t)pa[2] * pa[1] + (int128_t)pa[3] * pa[0]) << 1) + (((int128_t)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128_t)pa[3] * pa[1] + (int128_t)pa[4] * pa[0]) << 1) + (int128_t)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and PUTS the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[5];
	int128_t tmpZero[5];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi), with result coeffs centered
	tmpQ[0] = (int64_t)op[0] * -0x4a9b1119b4296ce4L + (int64_t)op[1] * 0x61109f06cb738ecaL + (int64_t)op[2] * 0x20706a64d7a3bb1fL + (int64_t)op[3] * 0x44585c239d88ea9dL + (int64_t)op[4] * -0x5cf825333b5e0f38L;
	tmpQ[1] = (int64_t)op[0] * 0x459be41eee90d4c2L + (int64_t)op[1] * -0x72e0a09c5ca274b4L + (int64_t)op[2] * 0x1b01b5fc4dcd694fL + (int64_t)op[3] * 0x2acce8c63958fdb6L + (int64_t)op[4] * 0x4c32e28e8b4eadf1L;
	tmpQ[2] = (int64_t)op[0] * -0x7fb8a39d5d89e4d1L + (int64_t)op[1] * -0x7418c9d933b31958L + (int64_t)op[2] * 0x353aca2988c7e7a8L + (int64_t)op[3] * 0x778e578180ad80c2L + (int64_t)op[4] * 0x776e4c77bcb64633L;
	tmpQ[3] = (int64_t)op[0] * -0x3dc30e18b341fcc7L + (int64_t)op[1] * -0xfd4e123b30bee8fL + (int64_t)op[2] * 0x160289018834d12dL + (int64_t)op[3] * 0x18fce730a82a17ccL + (int64_t)op[4] * -0x785afd104f402da1L;
	tmpQ[4] = (int64_t)op[0] * 0x71a98deda60650aaL + (int64_t)op[1] * -0x3c8bfb69fca5d8bdL + (int64_t)op[2] * 0x11b56449fbcc80dfL + (int64_t)op[3] * -0x267fc40eacf0ad68L + (int64_t)op[4] * -0xcb6374f3633e8bdL;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * -0x2453387d6afd0L + (int128_t)tmpQ[1] * -0xea11d6374798L + (int128_t)tmpQ[2] * -0x4bdc85c8f3306L + (int128_t)tmpQ[3] * 0x51f0e83f52cdbL + (int128_t)tmpQ[4] * -0x78d86be19fe69L;
	tmpZero[1] = (int128_t)tmpQ[0] * 0x3324934d37cdaL + (int128_t)tmpQ[1] * 0x56ab6df6d8b51L + (int128_t)tmpQ[2] * 0xea11d6374798L + (int128_t)tmpQ[3] * 0x2e6a0d95b1e64L + (int128_t)tmpQ[4] * -0x142800eec4767L;
	tmpZero[2] = (int128_t)tmpQ[0] * 0xa3674e968dd8L + (int128_t)tmpQ[1] * 0x2df13a2257335L + (int128_t)tmpQ[2] * -0x56ab6df6d8b51L + (int128_t)tmpQ[3] * 0xaaf485cc3907L + (int128_t)tmpQ[4] * 0xae0399f0598b7L;
	tmpZero[3] = (int128_t)tmpQ[0] * 0x132a750f130ffL + (int128_t)tmpQ[1] * -0x496200b4dd5a9L + (int128_t)tmpQ[2] * -0x2df13a2257335L + (int128_t)tmpQ[3] * 0x673bbe20b3fafL + (int128_t)tmpQ[4] * 0x2fac825970915L;
	tmpZero[4] = (int128_t)tmpQ[0] * -0x43556f80417aaL + (int128_t)tmpQ[1] * 0xf2c1ac1ca3ceL + (int128_t)tmpQ[2] * 0x496200b4dd5a9L + (int128_t)tmpQ[3] * 0x14ba33dea8437L + (int128_t)tmpQ[4] * 0x7459a1fb2b2afL;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
}

//~ --------------------------------------------------------------------

void exact_coeffs_reduction(int64_t *rop, int64_t *op){

	int i;
	int128_t tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++){
		tmp[i] = (int128_t) op[i];
	}

	internal_reduction(rop, tmp);

	mult_mod_poly(rop, rop, polys_P[0]);
}

//~ computes : op/phi
void from_mont_domain(int64_t *rop, int64_t *op){

	int i;
	int128_t tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++){
		tmp[i] = (int128_t) op[i];
	}

	internal_reduction(rop, tmp);
}

//~ return 0 if pA and pB represent the same element in Zp, and something else otherwise 
int64_t equality_check_polys(int64_t *pA, int64_t *pB){

	int i;
	int64_t res;
	int64_t poly_diff[NB_COEFF];
	int128_t tmp_result[NB_COEFF];

	sub_poly(poly_diff, pA, pB);

	for(i=0; i<NB_COEFF; i++){
		tmp_result[i] = (int128_t) poly_diff[i];
	}

	internal_reduction(poly_diff, tmp_result);

	res = 0;
	for(i=0; i<NB_COEFF; i++)
		res &= poly_diff[i];

	return res;
}

