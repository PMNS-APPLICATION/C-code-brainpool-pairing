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
//~ Note: E(X) = X^8 - X - 1
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];
	int128_t c8, c9, c10, c11, c12, c13, c14;

	c8 = (int128_t)pa[1] * pb[7] + (int128_t)pa[2] * pb[6] + (int128_t)pa[3] * pb[5] + (int128_t)pa[4] * pb[4] + (int128_t)pa[5] * pb[3] + (int128_t)pa[6] * pb[2] + (int128_t)pa[7] * pb[1];
	c9 = (int128_t)pa[2] * pb[7] + (int128_t)pa[3] * pb[6] + (int128_t)pa[4] * pb[5] + (int128_t)pa[5] * pb[4] + (int128_t)pa[6] * pb[3] + (int128_t)pa[7] * pb[2];
	c10 = (int128_t)pa[3] * pb[7] + (int128_t)pa[4] * pb[6] + (int128_t)pa[5] * pb[5] + (int128_t)pa[6] * pb[4] + (int128_t)pa[7] * pb[3];
	c11 = (int128_t)pa[4] * pb[7] + (int128_t)pa[5] * pb[6] + (int128_t)pa[6] * pb[5] + (int128_t)pa[7] * pb[4];
	c12 = (int128_t)pa[5] * pb[7] + (int128_t)pa[6] * pb[6] + (int128_t)pa[7] * pb[5];
	c13 = (int128_t)pa[6] * pb[7] + (int128_t)pa[7] * pb[6];
	c14 = (int128_t)pa[7] * pb[7];

	tmp_prod_result[0] = (int128_t)pa[0] * pb[0] + c8;
	tmp_prod_result[1] = (int128_t)pa[0] * pb[1] + (int128_t)pa[1] * pb[0] + c8 + c9;
	tmp_prod_result[2] = (int128_t)pa[0] * pb[2] + (int128_t)pa[1] * pb[1] + (int128_t)pa[2] * pb[0] + c9 + c10;
	tmp_prod_result[3] = (int128_t)pa[0] * pb[3] + (int128_t)pa[1] * pb[2] + (int128_t)pa[2] * pb[1] + (int128_t)pa[3] * pb[0] + c10 + c11;
	tmp_prod_result[4] = (int128_t)pa[0] * pb[4] + (int128_t)pa[1] * pb[3] + (int128_t)pa[2] * pb[2] + (int128_t)pa[3] * pb[1] + (int128_t)pa[4] * pb[0] + c11 + c12;
	tmp_prod_result[5] = (int128_t)pa[0] * pb[5] + (int128_t)pa[1] * pb[4] + (int128_t)pa[2] * pb[3] + (int128_t)pa[3] * pb[2] + (int128_t)pa[4] * pb[1] + (int128_t)pa[5] * pb[0] + c12 + c13;
	tmp_prod_result[6] = (int128_t)pa[0] * pb[6] + (int128_t)pa[1] * pb[5] + (int128_t)pa[2] * pb[4] + (int128_t)pa[3] * pb[3] + (int128_t)pa[4] * pb[2] + (int128_t)pa[5] * pb[1] + (int128_t)pa[6] * pb[0] + c13 + c14;
	tmp_prod_result[7] = (int128_t)pa[0] * pb[7] + (int128_t)pa[1] * pb[6] + (int128_t)pa[2] * pb[5] + (int128_t)pa[3] * pb[4] + (int128_t)pa[4] * pb[3] + (int128_t)pa[5] * pb[2] + (int128_t)pa[6] * pb[1] + (int128_t)pa[7] * pb[0] + c14;

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes: pa^2 mod(E)
//~ Note: E(X) = X^8 - X - 1
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128_t tmp_prod_result[NB_COEFF];
	int128_t c8, c9, c10, c11, c12, c13, c14;

	c8 = (((int128_t)pa[5] * pa[3] + (int128_t)pa[6] * pa[2] + (int128_t)pa[7] * pa[1]) << 1) + (int128_t)pa[4] * pa[4];
	c9 = (((int128_t)pa[5] * pa[4] + (int128_t)pa[6] * pa[3] + (int128_t)pa[7] * pa[2]) << 1);
	c10 = (((int128_t)pa[6] * pa[4] + (int128_t)pa[7] * pa[3]) << 1) + (int128_t)pa[5] * pa[5];
	c11 = (((int128_t)pa[6] * pa[5] + (int128_t)pa[7] * pa[4]) << 1);
	c12 = (((int128_t)pa[7] * pa[5]) << 1) + (int128_t)pa[6] * pa[6];
	c13 = (((int128_t)pa[7] * pa[6]) << 1);
	c14 = (int128_t)pa[7] * pa[7];

	tmp_prod_result[0] = (int128_t)pa[0] * pa[0] + c8;
	tmp_prod_result[1] = (((int128_t)pa[1] * pa[0]) << 1) + c8 + c9;
	tmp_prod_result[2] = (((int128_t)pa[2] * pa[0]) << 1) + (int128_t)pa[1] * pa[1] + c9 + c10;
	tmp_prod_result[3] = (((int128_t)pa[2] * pa[1] + (int128_t)pa[3] * pa[0]) << 1) + c10 + c11;
	tmp_prod_result[4] = (((int128_t)pa[3] * pa[1] + (int128_t)pa[4] * pa[0]) << 1) + (int128_t)pa[2] * pa[2] + c11 + c12;
	tmp_prod_result[5] = (((int128_t)pa[3] * pa[2] + (int128_t)pa[4] * pa[1] + (int128_t)pa[5] * pa[0]) << 1) + c12 + c13;
	tmp_prod_result[6] = (((int128_t)pa[4] * pa[2] + (int128_t)pa[5] * pa[1] + (int128_t)pa[6] * pa[0]) << 1) + (int128_t)pa[3] * pa[3] + c13 + c14;
	tmp_prod_result[7] = (((int128_t)pa[4] * pa[3] + (int128_t)pa[5] * pa[2] + (int128_t)pa[6] * pa[1] + (int128_t)pa[7] * pa[0]) << 1) + c14;

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and PUTS the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[8];
	int128_t tmpZero[8];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi), with result coeffs centered
	tmpQ[0] = (int64_t)op[0] * 0x11f73526aabcf22dL + (int64_t)op[1] * 0x16d3153ca2f72769L + (int64_t)op[2] * 0x56ea57e929d4a8d5L + (int64_t)op[3] * 0x7c71aa99dd220e39L + (int64_t)op[4] * -0x6597cddac8d829bdL + (int64_t)op[5] * 0x745ee5d7ffd61f7fL + (int64_t)op[6] * -0x7425280dd0cc0b24L + (int64_t)op[7] * 0x2f5dcf890675b45L;
	tmpQ[1] = (int64_t)op[0] * 0x186a55c1782f1f52L + (int64_t)op[1] * 0x4c94b4f9adf7b744L + (int64_t)op[2] * 0x72686f19b037378L + (int64_t)op[3] * -0x78ee399448a0f0dcL + (int64_t)op[4] * -0x2c73c045e9ce3edfL + (int64_t)op[5] * 0x508c23fa54da8a22L + (int64_t)op[6] * -0xe45453531ceab8bL + (int64_t)op[7] * -0x2de6d28d786749eaL;
	tmpQ[2] = (int64_t)op[0] * -0x47f28ac12b3e0b9fL + (int64_t)op[1] * -0x7f1bdb1fe81e91c5L + (int64_t)op[2] * -0x608d6ed8ec30cb3cL + (int64_t)op[3] * 0x2f6df95463e9e706L + (int64_t)op[4] * -0x579f7c8b2233cff6L + (int64_t)op[5] * 0x7383b84c132512f3L + (int64_t)op[6] * 0x5652732d87ecaa56L + (int64_t)op[7] * -0x399309e0aa4a5acdL;
	tmpQ[3] = (int64_t)op[0] * 0x6e8d34fef00aea8bL + (int64_t)op[1] * -0x32f17f0f420a09aL + (int64_t)op[2] * 0x63826d6d7bced1d9L + (int64_t)op[3] * 0x777b274cba07e5d1L + (int64_t)op[4] * 0x40e6a706bb1a8b53L + (int64_t)op[5] * 0x11745d2625ef4c25L + (int64_t)op[6] * -0x359b9ae0036787efL + (int64_t)op[7] * 0x2385cdf81fd14399L;
	tmpQ[4] = (int64_t)op[0] * 0x32f17f0f420a09aL + (int64_t)op[1] * -0x63826d6d7bced1d9L + (int64_t)op[2] * -0x777b274cba07e5d1L + (int64_t)op[3] * -0x40e6a706bb1a8b53L + (int64_t)op[4] * -0x11745d2625ef4c25L + (int64_t)op[5] * 0x359b9ae0036787efL + (int64_t)op[6] * -0x2385cdf81fd14399L + (int64_t)op[7] * -0x6b5e1d0dfbea49f1L;
	tmpQ[5] = (int64_t)op[0] * 0x7f1bdb1fe81e91c5L + (int64_t)op[1] * 0x608d6ed8ec30cb3cL + (int64_t)op[2] * -0x2f6df95463e9e706L + (int64_t)op[3] * 0x579f7c8b2233cff6L + (int64_t)op[4] * -0x7383b84c132512f3L + (int64_t)op[5] * -0x5652732d87ecaa56L + (int64_t)op[6] * 0x399309e0aa4a5acdL + (int64_t)op[7] * -0x38f19a1eeca3629cL;
	tmpQ[6] = (int64_t)op[0] * -0x3bf023b9980062ebL + (int64_t)op[1] * 0x480d2df8561dfecbL + (int64_t)op[2] * -0x6779dc6e22b1a4b7L + (int64_t)op[3] * -0x620f5b25ed35c6ceL + (int64_t)op[4] * 0x7411f1f274abcdbbL + (int64_t)op[5] * 0x5d18d7d8ca1b9e66L + (int64_t)op[6] * 0x326c82ef0f46e755L + (int64_t)op[7] * 0x3ffc9f755bfd8e40L;
	tmpQ[7] = (int64_t)op[0] * 0x328726263a26da81L + (int64_t)op[1] * 0x5966e7e7512d57c4L + (int64_t)op[2] * 0x4980403fe4b709d6L + (int64_t)op[3] * -0x7becc32ef3fdf12bL + (int64_t)op[4] * 0x3bf023b9980062ebL + (int64_t)op[5] * -0x480d2df8561dfecbL + (int64_t)op[6] * 0x6779dc6e22b1a4b7L + (int64_t)op[7] * 0x620f5b25ed35c6ceL;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * -0xb8cf6d636dea08L + (int128_t)tmpQ[1] * -0x2faef54891fa477L + (int128_t)tmpQ[2] * -0x2d4f32024b4bd92L + (int128_t)tmpQ[3] * 0x78ded52f33497cL + (int128_t)tmpQ[4] * -0x16853fcf36aa554L + (int128_t)tmpQ[5] * 0x612f3383c31507L + (int128_t)tmpQ[6] * -0x14e22bc4be39a5aL + (int128_t)tmpQ[7] * 0xdc760631b3aa0aL;
	tmpZero[1] = (int128_t)tmpQ[0] * 0x21c23b2c146d38aL + (int128_t)tmpQ[1] * -0xc281a261e3f063L + (int128_t)tmpQ[2] * -0x117a905510bd82cL + (int128_t)tmpQ[3] * -0x1929b5795b4ff23L + (int128_t)tmpQ[4] * -0x402d750bf1f1dfL + (int128_t)tmpQ[5] * -0xecf388c8208553L + (int128_t)tmpQ[6] * 0x1c6c0fc0d0bc3c4L + (int128_t)tmpQ[7] * 0x195457395219412L;
	tmpZero[2] = (int128_t)tmpQ[0] * 0x117a905510bd82cL + (int128_t)tmpQ[1] * -0x10724c96fa7904dL + (int128_t)tmpQ[2] * 0x1a2fb794e2e1a4L + (int128_t)tmpQ[3] * 0x110472a3fc3009fL + (int128_t)tmpQ[4] * -0x17199991de3dbfbL + (int128_t)tmpQ[5] * 0x1c6c0fc0d0bc3c4L + (int128_t)tmpQ[6] * 0x8e20aa257a03c5L + (int128_t)tmpQ[7] * -0x21c23b2c146d38aL;
	tmpZero[3] = (int128_t)tmpQ[0] * -0x1a2fb794e2e1a4L + (int128_t)tmpQ[1] * -0x25fc34646ae6e5L + (int128_t)tmpQ[2] * 0x15c878023e2a013L + (int128_t)tmpQ[3] * -0xa5d25a9186b4f1L + (int128_t)tmpQ[4] * -0x1af752b2a0475b6L + (int128_t)tmpQ[5] * 0x8e20aa257a03c5L + (int128_t)tmpQ[6] * -0x2421fe725b1ba6fL + (int128_t)tmpQ[7] * -0x117a905510bd82cL;
	tmpZero[4] = (int128_t)tmpQ[0] * -0x15c878023e2a013L + (int128_t)tmpQ[1] * 0x7b23975392ceaeL + (int128_t)tmpQ[2] * 0xaf848f8ffcbb4cL + (int128_t)tmpQ[3] * 0x22f4b51572043c2L + (int128_t)tmpQ[4] * 0xbf9c329d643377L + (int128_t)tmpQ[5] * -0x2421fe725b1ba6fL + (int128_t)tmpQ[6] * -0x9c856dfd79097eL + (int128_t)tmpQ[7] * 0x1a2fb794e2e1a4L;
	tmpZero[5] = (int128_t)tmpQ[0] * -0xaf848f8ffcbb4cL + (int128_t)tmpQ[1] * -0x19c7818581d409fL + (int128_t)tmpQ[2] * 0x1accc983d3c0a1dL + (int128_t)tmpQ[3] * -0x273c3eca0f1a88bL + (int128_t)tmpQ[4] * 0x682475c10f1ce0L + (int128_t)tmpQ[5] * -0x9c856dfd79097eL + (int128_t)tmpQ[6] * -0x1824860c33a5efbL + (int128_t)tmpQ[7] * 0x15c878023e2a013L;
	tmpZero[6] = (int128_t)tmpQ[0] * -0x1accc983d3c0a1dL + (int128_t)tmpQ[1] * 0x19f463cfcfb9a7L + (int128_t)tmpQ[2] * 0xdc760631b3aa0aL + (int128_t)tmpQ[3] * 0x6f275e87c54b0cL + (int128_t)tmpQ[4] * -0x1c6fc4fd21eebc1L + (int128_t)tmpQ[5] * -0x1824860c33a5efbL + (int128_t)tmpQ[6] * 0x1767be3f3b259baL + (int128_t)tmpQ[7] * 0xaf848f8ffcbb4cL;
	tmpZero[7] = (int128_t)tmpQ[0] * -0xdc760631b3aa0aL + (int128_t)tmpQ[1] * -0x4e555c0c39a645L + (int128_t)tmpQ[2] * 0xb8cf6d636dea08L + (int128_t)tmpQ[3] * 0x171c9551a295a5cL + (int128_t)tmpQ[4] * -0x238fd8655964a1dL + (int128_t)tmpQ[5] * 0x1767be3f3b259baL + (int128_t)tmpQ[6] * 0x612f3383c31507L + (int128_t)tmpQ[7] * 0x1accc983d3c0a1dL;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
	rop[5] = (op[5] + tmpZero[5]) >> PHI_LOG2;
	rop[6] = (op[6] + tmpZero[6]) >> PHI_LOG2;
	rop[7] = (op[7] + tmpZero[7]) >> PHI_LOG2;
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

