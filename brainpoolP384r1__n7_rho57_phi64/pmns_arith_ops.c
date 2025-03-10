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
//~ Note: E(X) = X^7 - 2
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pb[0] + (((int128_t)pa[1] * pb[6] + (int128_t)pa[2] * pb[5] + (int128_t)pa[3] * pb[4] + (int128_t)pa[4] * pb[3] + (int128_t)pa[5] * pb[2] + (int128_t)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128_t)pa[0] * pb[1] + (int128_t)pa[1] * pb[0] + (((int128_t)pa[2] * pb[6] + (int128_t)pa[3] * pb[5] + (int128_t)pa[4] * pb[4] + (int128_t)pa[5] * pb[3] + (int128_t)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128_t)pa[0] * pb[2] + (int128_t)pa[1] * pb[1] + (int128_t)pa[2] * pb[0] + (((int128_t)pa[3] * pb[6] + (int128_t)pa[4] * pb[5] + (int128_t)pa[5] * pb[4] + (int128_t)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128_t)pa[0] * pb[3] + (int128_t)pa[1] * pb[2] + (int128_t)pa[2] * pb[1] + (int128_t)pa[3] * pb[0] + (((int128_t)pa[4] * pb[6] + (int128_t)pa[5] * pb[5] + (int128_t)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128_t)pa[0] * pb[4] + (int128_t)pa[1] * pb[3] + (int128_t)pa[2] * pb[2] + (int128_t)pa[3] * pb[1] + (int128_t)pa[4] * pb[0] + (((int128_t)pa[5] * pb[6] + (int128_t)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128_t)pa[0] * pb[5] + (int128_t)pa[1] * pb[4] + (int128_t)pa[2] * pb[3] + (int128_t)pa[3] * pb[2] + (int128_t)pa[4] * pb[1] + (int128_t)pa[5] * pb[0] + (((int128_t)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128_t)pa[0] * pb[6] + (int128_t)pa[1] * pb[5] + (int128_t)pa[2] * pb[4] + (int128_t)pa[3] * pb[3] + (int128_t)pa[4] * pb[2] + (int128_t)pa[5] * pb[1] + (int128_t)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes: pa^2 mod(E)
//~ Note: E(X) = X^7 - 2
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pa[0] + (((int128_t)pa[4] * pa[3] + (int128_t)pa[5] * pa[2] + (int128_t)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128_t)pa[1] * pa[0]) << 1) + (((((int128_t)pa[5] * pa[3] + (int128_t)pa[6] * pa[2]) << 1) + (int128_t)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128_t)pa[2] * pa[0]) << 1) + (int128_t)pa[1] * pa[1] + (((int128_t)pa[5] * pa[4] + (int128_t)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128_t)pa[2] * pa[1] + (int128_t)pa[3] * pa[0]) << 1) + (((((int128_t)pa[6] * pa[4]) << 1) + (int128_t)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128_t)pa[3] * pa[1] + (int128_t)pa[4] * pa[0]) << 1) + (int128_t)pa[2] * pa[2] + (((int128_t)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128_t)pa[3] * pa[2] + (int128_t)pa[4] * pa[1] + (int128_t)pa[5] * pa[0]) << 1) + (((int128_t)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128_t)pa[4] * pa[2] + (int128_t)pa[5] * pa[1] + (int128_t)pa[6] * pa[0]) << 1) + (int128_t)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and PUTS the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[7];
	int128_t tmpZero[7];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi), with result coeffs centered
	tmpQ[0] = (int64_t)op[0] * 0x4576b194143be862L + (int64_t)op[1] * -0x4117eb301436a058L + (int64_t)op[2] * 0x15400691f655124aL + (int64_t)op[3] * 0x485b505c3e7a035dL + (int64_t)op[4] * -0x7de17dbfe1de96bfL + (int64_t)op[5] * 0x5dbcef3a741e76f0L + (int64_t)op[6] * 0x66a1327424d19e2aL;
	tmpQ[1] = (int64_t)op[0] * -0x5ad9cdfaa6928487L + (int64_t)op[1] * -0x912e08e5f54f551L + (int64_t)op[2] * -0x66f2e59d6f74f61cL + (int64_t)op[3] * 0x309c6c8c36a7d96eL + (int64_t)op[4] * -0x29e553674cce71e7L + (int64_t)op[5] * 0x2eff2f1bc4df2284L + (int64_t)op[6] * -0x3ca67cbb31a66c8cL;
	tmpQ[2] = (int64_t)op[0] * -0x57f52f327c54698cL + (int64_t)op[1] * -0x45dc731e2cfcebb8L + (int64_t)op[2] * -0x1e75fcf4f1ab9176L + (int64_t)op[3] * 0x4ee24ea41cff743bL + (int64_t)op[4] * -0x2116727f42780a64L + (int64_t)op[5] * 0x4f12698128536ae4L + (int64_t)op[6] * -0x78c7a20b69cde622L;
	tmpQ[3] = (int64_t)op[0] * -0x45dc731e2cfcebb8L + (int64_t)op[1] * -0x1e75fcf4f1ab9176L + (int64_t)op[2] * 0x4ee24ea41cff743bL + (int64_t)op[3] * -0x2116727f42780a64L + (int64_t)op[4] * 0x4f12698128536ae4L + (int64_t)op[5] * -0x78c7a20b69cde622L + (int64_t)op[6] * 0x5015a19b07572ce8L;
	tmpQ[4] = (int64_t)op[0] * 0x3c63d105b4e6f311L + (int64_t)op[1] * 0x57f52f327c54698cL + (int64_t)op[2] * 0x45dc731e2cfcebb8L + (int64_t)op[3] * 0x1e75fcf4f1ab9176L + (int64_t)op[4] * -0x4ee24ea41cff743bL + (int64_t)op[5] * 0x2116727f42780a64L + (int64_t)op[6] * -0x4f12698128536ae4L;
	tmpQ[5] = (int64_t)op[0] * -0x278934c09429b572L + (int64_t)op[1] * 0x3c63d105b4e6f311L + (int64_t)op[2] * 0x57f52f327c54698cL + (int64_t)op[3] * 0x45dc731e2cfcebb8L + (int64_t)op[4] * 0x1e75fcf4f1ab9176L + (int64_t)op[5] * -0x4ee24ea41cff743bL + (int64_t)op[6] * 0x2116727f42780a64L;
	tmpQ[6] = (int64_t)op[0] * -0x108b393fa13c0532L + (int64_t)op[1] * 0x278934c09429b572L + (int64_t)op[2] * -0x3c63d105b4e6f311L + (int64_t)op[3] * -0x57f52f327c54698cL + (int64_t)op[4] * -0x45dc731e2cfcebb8L + (int64_t)op[5] * -0x1e75fcf4f1ab9176L + (int64_t)op[6] * 0x4ee24ea41cff743bL;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * 0x5783b6a3a90fabL + (int128_t)tmpQ[1] * 0x226dc9ae6a0484L + (int128_t)tmpQ[2] * 0x1ae6f50cdcf1a9L + (int128_t)tmpQ[3] * 0x392e4ef22ce648L + (int128_t)tmpQ[4] * -0x283438ea027fL + (int128_t)tmpQ[5] * 0x17210a5c10b798L + (int128_t)tmpQ[6] * 0x791a76fd96308L;
	tmpZero[1] = (int128_t)tmpQ[0] * -0x3c8d3b7ecb184L + (int128_t)tmpQ[1] * -0x5783b6a3a90fabL + (int128_t)tmpQ[2] * 0x35657b3a4034c4L + (int128_t)tmpQ[3] * 0x226dc9ae6a0484L + (int128_t)tmpQ[4] * -0x171e2154f04025L + (int128_t)tmpQ[5] * 0x34edb8bc5508a8L + (int128_t)tmpQ[6] * -0x17210a5c10b798L;
	tmpZero[2] = (int128_t)tmpQ[0] * 0xb90852e085bccL + (int128_t)tmpQ[1] * 0x3c8d3b7ecb184L + (int128_t)tmpQ[2] * 0x2dfe4edc726050L + (int128_t)tmpQ[3] * -0x5783b6a3a90fabL + (int128_t)tmpQ[4] * -0x40f60068489090L + (int128_t)tmpQ[5] * 0x3c9cc196cc1e02L + (int128_t)tmpQ[6] * -0x34edb8bc5508a8L;
	tmpZero[3] = (int128_t)tmpQ[0] * 0x1a76dc5e2a8454L + (int128_t)tmpQ[1] * -0xb90852e085bccL + (int128_t)tmpQ[2] * -0x3d0cda457e8b57L + (int128_t)tmpQ[3] * 0x3c8d3b7ecb184L + (int128_t)tmpQ[4] * -0x48752b3a9ce4a4L + (int128_t)tmpQ[5] * -0x392e4ef22ce648L + (int128_t)tmpQ[6] * -0x3c9cc196cc1e02L;
	tmpZero[4] = (int128_t)tmpQ[0] * 0x1e4e60cb660f01L + (int128_t)tmpQ[1] * -0x1a76dc5e2a8454L + (int128_t)tmpQ[2] * 0x2217348352c085L + (int128_t)tmpQ[3] * -0xb90852e085bccL + (int128_t)tmpQ[4] * 0x1ebe797a187c56L + (int128_t)tmpQ[5] * -0x226dc9ae6a0484L + (int128_t)tmpQ[6] * 0x392e4ef22ce648L;
	tmpZero[5] = (int128_t)tmpQ[0] * -0x1c972779167324L + (int128_t)tmpQ[1] * -0x1e4e60cb660f01L + (int128_t)tmpQ[2] * -0x2827aca71ecef0L + (int128_t)tmpQ[3] * -0x1a76dc5e2a8454L + (int128_t)tmpQ[4] * -0x5800d0a3c4d61L + (int128_t)tmpQ[5] * 0x5783b6a3a90fabL + (int128_t)tmpQ[6] * 0x226dc9ae6a0484L;
	tmpZero[6] = (int128_t)tmpQ[0] * -0x1136e4d7350242L + (int128_t)tmpQ[1] * 0x1c972779167324L + (int128_t)tmpQ[2] * -0x2badc1355f8696L + (int128_t)tmpQ[3] * -0x1e4e60cb660f01L + (int128_t)tmpQ[4] * 0x395e917e53d132L + (int128_t)tmpQ[5] * -0x3c8d3b7ecb184L + (int128_t)tmpQ[6] * -0x5783b6a3a90fabL;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
	rop[5] = (op[5] + tmpZero[5]) >> PHI_LOG2;
	rop[6] = (op[6] + tmpZero[6]) >> PHI_LOG2;
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

