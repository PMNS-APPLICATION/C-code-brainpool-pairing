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
	tmpQ[0] = (int64_t)op[0] * 0x19472d32a6c893cL + (int64_t)op[1] * -0x7f51c29809ee2999L + (int64_t)op[2] * -0x267151b975af6de8L + (int64_t)op[3] * -0xc1614faa2d035afL + (int64_t)op[4] * 0x5d92bfe52c20abffL + (int64_t)op[5] * -0x56cae0cc12987c2fL + (int64_t)op[6] * -0x7ebd0fe79a18619dL;
	tmpQ[1] = (int64_t)op[0] * -0x7657f924130ed0ceL + (int64_t)op[1] * 0xa1e69575fdc9426L + (int64_t)op[2] * -0x1f5fa48f95b86754L + (int64_t)op[3] * 0x2e6b3df3f6d59a78L + (int64_t)op[4] * 0x1b0da555eb47bd11L + (int64_t)op[5] * 0x776bc24fe64b3dd0L + (int64_t)op[6] * -0x3e90023a64a1d065L;
	tmpQ[2] = (int64_t)op[0] * -0x67fd0d3227c95d5aL + (int64_t)op[1] * -0x15c97e6a74617c9L + (int64_t)op[2] * -0x7da1813a762c48d5L + (int64_t)op[3] * -0x336fa5c4b9aec993L + (int64_t)op[4] * 0x7d066d7dda7ebf15L + (int64_t)op[5] * -0x10bbbdab39bca8edL + (int64_t)op[6] * -0x71f66a913c92d739L;
	tmpQ[3] = (int64_t)op[0] * -0x7840156e29f72305L + (int64_t)op[1] * 0x1410013523f6623fL + (int64_t)op[2] * -0x4e9b2f89e3a9249dL + (int64_t)op[3] * 0x2bc96301250465feL + (int64_t)op[4] * -0x169dd31edd21eaf7L + (int64_t)op[5] * -0x6e95e7d6150f15b1L + (int64_t)op[6] * 0x16093d85286e8df6L;
	tmpQ[4] = (int64_t)op[0] * 0x13acc8650b6f49f2L + (int64_t)op[1] * 0x70462933a8d0f675L + (int64_t)op[2] * -0x230791906c6de9a1L + (int64_t)op[3] * 0x799cad54d13d3aceL + (int64_t)op[4] * 0x699180055235e696L + (int64_t)op[5] * -0x7df52ab162a811d0L + (int64_t)op[6] * 0x57302f81007cdaedL;
	tmpQ[5] = (int64_t)op[0] * 0x74fb613d6bc8b905L + (int64_t)op[1] * 0x7840156e29f72305L + (int64_t)op[2] * -0x1410013523f6623fL + (int64_t)op[3] * 0x4e9b2f89e3a9249dL + (int64_t)op[4] * -0x2bc96301250465feL + (int64_t)op[5] * 0x169dd31edd21eaf7L + (int64_t)op[6] * 0x6e95e7d6150f15b1L;
	tmpQ[6] = (int64_t)op[0] * -0xbb2dc2a8a491d62L + (int64_t)op[1] * -0x614e98d860596f13L + (int64_t)op[2] * -0x7f9ec3a81262c90L + (int64_t)op[3] * -0xef7905b48778762L + (int64_t)op[4] * 0x2b017dcaed941631L + (int64_t)op[5] * -0x6aa51cf988c5b36cL + (int64_t)op[6] * 0x6b6d022fc0360339L;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * 0xf56519856382dL + (int128_t)tmpQ[1] * -0x1332c7b9665d06L + (int128_t)tmpQ[2] * -0xaf5b3ddf9c5eL + (int128_t)tmpQ[3] * 0x1ab484f69bc266L + (int128_t)tmpQ[4] * 0x12e343a51b769L + (int128_t)tmpQ[5] * 0x45d8697a692ba1L + (int128_t)tmpQ[6] * -0x3b0662b5c57025L;
	tmpZero[1] = (int128_t)tmpQ[0] * 0x51f35f3f2ea602L + (int128_t)tmpQ[1] * -0xf56519856382dL + (int128_t)tmpQ[2] * -0x1332c7b9665d06L + (int128_t)tmpQ[3] * -0xaf5b3ddf9c5eL + (int128_t)tmpQ[4] * -0x846062ef29a0eL + (int128_t)tmpQ[5] * -0x131dda0588c4b6L + (int128_t)tmpQ[6] * -0x717d1f4a0e2a5L;
	tmpZero[2] = (int128_t)tmpQ[0] * 0x1ce56d6f220285L + (int128_t)tmpQ[1] * -0x51f35f3f2ea602L + (int128_t)tmpQ[2] * -0xf56519856382dL + (int128_t)tmpQ[3] * -0x1332c7b9665d06L + (int128_t)tmpQ[4] * 0x28f2fc0b47291cL + (int128_t)tmpQ[5] * -0x3b766886cde9c4L + (int128_t)tmpQ[6] * 0x20acf5dc548f0eL;
	tmpZero[3] = (int128_t)tmpQ[0] * 0x16290c5f0f9052L + (int128_t)tmpQ[1] * -0x1ce56d6f220285L + (int128_t)tmpQ[2] * -0x51f35f3f2ea602L + (int128_t)tmpQ[3] * -0xf56519856382dL + (int128_t)tmpQ[4] * -0x2946e664985508L + (int128_t)tmpQ[5] * 0x2d235c85a879e1L + (int128_t)tmpQ[6] * -0x53ea59512becL;
	tmpZero[4] = (int128_t)tmpQ[0] * -0xd5a427b4de133L + (int128_t)tmpQ[1] * -0x16290c5f0f9052L + (int128_t)tmpQ[2] * -0x1ce56d6f220285L + (int128_t)tmpQ[3] * -0x51f35f3f2ea602L + (int128_t)tmpQ[4] * -0x2e1c260b800891L + (int128_t)tmpQ[5] * -0x1480e79b586544L + (int128_t)tmpQ[6] * -0x57630c70185d99L;
	tmpZero[5] = (int128_t)tmpQ[0] * 0x57ad9eefce2fL + (int128_t)tmpQ[1] * 0xd5a427b4de133L + (int128_t)tmpQ[2] * -0x16290c5f0f9052L + (int128_t)tmpQ[3] * -0x1ce56d6f220285L + (int128_t)tmpQ[4] * 0x2ccbaee6b8abb2L + (int128_t)tmpQ[5] * 0x84242e953f7cbL + (int128_t)tmpQ[6] * -0x1507724c75cdfL;
	tmpZero[6] = (int128_t)tmpQ[0] * 0x99963dcb32e83L + (int128_t)tmpQ[1] * -0x57ad9eefce2fL + (int128_t)tmpQ[2] * 0xd5a427b4de133L + (int128_t)tmpQ[3] * -0x16290c5f0f9052L + (int128_t)tmpQ[4] * -0x1e1a4b780b93c7L + (int128_t)tmpQ[5] * 0x24d6ac881e05faL + (int128_t)tmpQ[6] * 0xeb1636ead17ebL;

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

