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
//~ Note: E(X) = X^6 - X^3 + 1
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];
	int128_t c6, c7, c8, c9, c10;

	c6 = (int128_t)pa[1] * pb[5] + (int128_t)pa[2] * pb[4] + (int128_t)pa[3] * pb[3] + (int128_t)pa[4] * pb[2] + (int128_t)pa[5] * pb[1];
	c7 = (int128_t)pa[2] * pb[5] + (int128_t)pa[3] * pb[4] + (int128_t)pa[4] * pb[3] + (int128_t)pa[5] * pb[2];
	c8 = (int128_t)pa[3] * pb[5] + (int128_t)pa[4] * pb[4] + (int128_t)pa[5] * pb[3];
	c9 = (int128_t)pa[4] * pb[5] + (int128_t)pa[5] * pb[4];
	c10 = (int128_t)pa[5] * pb[5];

	tmp_prod_result[0] = (int128_t)pa[0] * pb[0] - c6 - c9;
	tmp_prod_result[1] = (int128_t)pa[0] * pb[1] + (int128_t)pa[1] * pb[0] - c7 - c10;
	tmp_prod_result[2] = (int128_t)pa[0] * pb[2] + (int128_t)pa[1] * pb[1] + (int128_t)pa[2] * pb[0] - c8;
	tmp_prod_result[3] = (int128_t)pa[0] * pb[3] + (int128_t)pa[1] * pb[2] + (int128_t)pa[2] * pb[1] + (int128_t)pa[3] * pb[0] + c6;
	tmp_prod_result[4] = (int128_t)pa[0] * pb[4] + (int128_t)pa[1] * pb[3] + (int128_t)pa[2] * pb[2] + (int128_t)pa[3] * pb[1] + (int128_t)pa[4] * pb[0] + c7;
	tmp_prod_result[5] = (int128_t)pa[0] * pb[5] + (int128_t)pa[1] * pb[4] + (int128_t)pa[2] * pb[3] + (int128_t)pa[3] * pb[2] + (int128_t)pa[4] * pb[1] + (int128_t)pa[5] * pb[0] + c8;

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes: pa^2 mod(E)
//~ Note: E(X) = X^6 - X^3 + 1
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128_t tmp_prod_result[NB_COEFF];
	int128_t c6, c7, c8, c9, c10;

	c6 = (((int128_t)pa[4] * pa[2] + (int128_t)pa[5] * pa[1]) << 1) + (int128_t)pa[3] * pa[3];
	c7 = (((int128_t)pa[4] * pa[3] + (int128_t)pa[5] * pa[2]) << 1);
	c8 = (((int128_t)pa[5] * pa[3]) << 1) + (int128_t)pa[4] * pa[4];
	c9 = (((int128_t)pa[5] * pa[4]) << 1);
	c10 = (int128_t)pa[5] * pa[5];

	tmp_prod_result[0] = (int128_t)pa[0] * pa[0] - c6 - c9;
	tmp_prod_result[1] = (((int128_t)pa[1] * pa[0]) << 1) - c7 - c10;
	tmp_prod_result[2] = (((int128_t)pa[2] * pa[0]) << 1) + (int128_t)pa[1] * pa[1] - c8;
	tmp_prod_result[3] = (((int128_t)pa[2] * pa[1] + (int128_t)pa[3] * pa[0]) << 1) + c6;
	tmp_prod_result[4] = (((int128_t)pa[3] * pa[1] + (int128_t)pa[4] * pa[0]) << 1) + (int128_t)pa[2] * pa[2] + c7;
	tmp_prod_result[5] = (((int128_t)pa[3] * pa[2] + (int128_t)pa[4] * pa[1] + (int128_t)pa[5] * pa[0]) << 1) + c8;

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and PUTS the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[6];
	int128_t tmpZero[6];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi), with result coeffs centered
	tmpQ[0] = (int64_t)op[0] * -0x47c7ffa99bac3ca7L + (int64_t)op[1] * -0x4f999a9848f23609L + (int64_t)op[2] * -0x349d599a6be51eb8L + (int64_t)op[3] * -0x5c91d61b917c6bddL + (int64_t)op[4] * 0x75a364550899f9e6L + (int64_t)op[5] * -0x540adba9b7cf2afbL;
	tmpQ[1] = (int64_t)op[0] * -0x4f999a9848f23609L + (int64_t)op[1] * -0x349d599a6be51eb8L + (int64_t)op[2] * -0x5c91d61b917c6bddL + (int64_t)op[3] * 0x75a364550899f9e6L + (int64_t)op[4] * -0x540adba9b7cf2afbL + (int64_t)op[5] * -0x14c9d671f5d02f36L;
	tmpQ[2] = (int64_t)op[0] * 0x349d599a6be51eb8L + (int64_t)op[1] * 0x5c91d61b917c6bddL + (int64_t)op[2] * -0x75a364550899f9e6L + (int64_t)op[3] * 0x540adba9b7cf2afbL + (int64_t)op[4] * 0x14c9d671f5d02f36L + (int64_t)op[5] * 0x3ac30112ae73d011L;
	tmpQ[3] = (int64_t)op[0] * 0x75a364550899f9e6L + (int64_t)op[1] * -0x540adba9b7cf2afbL + (int64_t)op[2] * -0x14c9d671f5d02f36L + (int64_t)op[3] * -0x3ac30112ae73d011L + (int64_t)op[4] * -0x1f6d820f4bea0c43L + (int64_t)op[5] * 0x47c7ffa99bac3ca7L;
	tmpQ[4] = (int64_t)op[0] * 0x540adba9b7cf2afbL + (int64_t)op[1] * 0x14c9d671f5d02f36L + (int64_t)op[2] * 0x3ac30112ae73d011L + (int64_t)op[3] * 0x1f6d820f4bea0c43L + (int64_t)op[4] * -0x47c7ffa99bac3ca7L + (int64_t)op[5] * -0x4f999a9848f23609L;
	tmpQ[5] = (int64_t)op[0] * 0x5c91d61b917c6bddL + (int64_t)op[1] * -0x75a364550899f9e6L + (int64_t)op[2] * 0x540adba9b7cf2afbL + (int64_t)op[3] * 0x14c9d671f5d02f36L + (int64_t)op[4] * 0x3ac30112ae73d011L + (int64_t)op[5] * 0x1f6d820f4bea0c43L;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * -0x21edce7307bd6fL + (int128_t)tmpQ[1] * -0x29e6f55707a2cL + (int128_t)tmpQ[2] * -0x6201304c78be89L + (int128_t)tmpQ[3] * 0x14fe89a27bfd8aL + (int128_t)tmpQ[4] * -0x138e797826dd12L + (int128_t)tmpQ[5] * 0x10eb06fe05021cL;
	tmpZero[1] = (int128_t)tmpQ[0] * -0x29e6f55707a2cL + (int128_t)tmpQ[1] * 0x6201304c78be89L + (int128_t)tmpQ[2] * 0x32d8d5710cbf8bL + (int128_t)tmpQ[3] * 0x138e797826dd12L + (int128_t)tmpQ[4] * -0x21edce7307bd6fL + (int128_t)tmpQ[5] * -0x14fe89a27bfd8aL;
	tmpZero[2] = (int128_t)tmpQ[0] * 0x6201304c78be89L + (int128_t)tmpQ[1] * -0x32d8d5710cbf8bL + (int128_t)tmpQ[2] * -0x12601a4d0b835eL + (int128_t)tmpQ[3] * 0x21edce7307bd6fL + (int128_t)tmpQ[4] * -0x29e6f55707a2cL + (int128_t)tmpQ[5] * -0x138e797826dd12L;
	tmpZero[3] = (int128_t)tmpQ[0] * -0x10eb06fe05021cL + (int128_t)tmpQ[1] * 0x14fe89a27bfd8aL + (int128_t)tmpQ[2] * -0x138e797826dd12L + (int128_t)tmpQ[3] * -0x12601a4d0b835eL + (int128_t)tmpQ[4] * 0x758fa9c49f9b9bL + (int128_t)tmpQ[5] * -0x32d8d5710cbf8bL;
	tmpZero[4] = (int128_t)tmpQ[0] * 0x14fe89a27bfd8aL + (int128_t)tmpQ[1] * 0x138e797826dd12L + (int128_t)tmpQ[2] * -0x21edce7307bd6fL + (int128_t)tmpQ[3] * -0x758fa9c49f9b9bL + (int128_t)tmpQ[4] * -0x10eb06fe05021cL + (int128_t)tmpQ[5] * 0x12601a4d0b835eL;
	tmpZero[5] = (int128_t)tmpQ[0] * 0x138e797826dd12L + (int128_t)tmpQ[1] * 0x21edce7307bd6fL + (int128_t)tmpQ[2] * -0x29e6f55707a2cL + (int128_t)tmpQ[3] * 0x10eb06fe05021cL + (int128_t)tmpQ[4] * 0x14fe89a27bfd8aL + (int128_t)tmpQ[5] * 0x758fa9c49f9b9bL;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
	rop[5] = (op[5] + tmpZero[5]) >> PHI_LOG2;
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

