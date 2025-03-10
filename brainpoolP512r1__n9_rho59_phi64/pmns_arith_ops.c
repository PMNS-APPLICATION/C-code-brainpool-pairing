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
//~ Note: E(X) = X^9 - 2
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pb[0] + (((int128_t)pa[1] * pb[8] + (int128_t)pa[2] * pb[7] + (int128_t)pa[3] * pb[6] + (int128_t)pa[4] * pb[5] + (int128_t)pa[5] * pb[4] + (int128_t)pa[6] * pb[3] + (int128_t)pa[7] * pb[2] + (int128_t)pa[8] * pb[1]) << 1);
	tmp_prod_result[1] = (int128_t)pa[0] * pb[1] + (int128_t)pa[1] * pb[0] + (((int128_t)pa[2] * pb[8] + (int128_t)pa[3] * pb[7] + (int128_t)pa[4] * pb[6] + (int128_t)pa[5] * pb[5] + (int128_t)pa[6] * pb[4] + (int128_t)pa[7] * pb[3] + (int128_t)pa[8] * pb[2]) << 1);
	tmp_prod_result[2] = (int128_t)pa[0] * pb[2] + (int128_t)pa[1] * pb[1] + (int128_t)pa[2] * pb[0] + (((int128_t)pa[3] * pb[8] + (int128_t)pa[4] * pb[7] + (int128_t)pa[5] * pb[6] + (int128_t)pa[6] * pb[5] + (int128_t)pa[7] * pb[4] + (int128_t)pa[8] * pb[3]) << 1);
	tmp_prod_result[3] = (int128_t)pa[0] * pb[3] + (int128_t)pa[1] * pb[2] + (int128_t)pa[2] * pb[1] + (int128_t)pa[3] * pb[0] + (((int128_t)pa[4] * pb[8] + (int128_t)pa[5] * pb[7] + (int128_t)pa[6] * pb[6] + (int128_t)pa[7] * pb[5] + (int128_t)pa[8] * pb[4]) << 1);
	tmp_prod_result[4] = (int128_t)pa[0] * pb[4] + (int128_t)pa[1] * pb[3] + (int128_t)pa[2] * pb[2] + (int128_t)pa[3] * pb[1] + (int128_t)pa[4] * pb[0] + (((int128_t)pa[5] * pb[8] + (int128_t)pa[6] * pb[7] + (int128_t)pa[7] * pb[6] + (int128_t)pa[8] * pb[5]) << 1);
	tmp_prod_result[5] = (int128_t)pa[0] * pb[5] + (int128_t)pa[1] * pb[4] + (int128_t)pa[2] * pb[3] + (int128_t)pa[3] * pb[2] + (int128_t)pa[4] * pb[1] + (int128_t)pa[5] * pb[0] + (((int128_t)pa[6] * pb[8] + (int128_t)pa[7] * pb[7] + (int128_t)pa[8] * pb[6]) << 1);
	tmp_prod_result[6] = (int128_t)pa[0] * pb[6] + (int128_t)pa[1] * pb[5] + (int128_t)pa[2] * pb[4] + (int128_t)pa[3] * pb[3] + (int128_t)pa[4] * pb[2] + (int128_t)pa[5] * pb[1] + (int128_t)pa[6] * pb[0] + (((int128_t)pa[7] * pb[8] + (int128_t)pa[8] * pb[7]) << 1);
	tmp_prod_result[7] = (int128_t)pa[0] * pb[7] + (int128_t)pa[1] * pb[6] + (int128_t)pa[2] * pb[5] + (int128_t)pa[3] * pb[4] + (int128_t)pa[4] * pb[3] + (int128_t)pa[5] * pb[2] + (int128_t)pa[6] * pb[1] + (int128_t)pa[7] * pb[0] + (((int128_t)pa[8] * pb[8]) << 1);
	tmp_prod_result[8] = (int128_t)pa[0] * pb[8] + (int128_t)pa[1] * pb[7] + (int128_t)pa[2] * pb[6] + (int128_t)pa[3] * pb[5] + (int128_t)pa[4] * pb[4] + (int128_t)pa[5] * pb[3] + (int128_t)pa[6] * pb[2] + (int128_t)pa[7] * pb[1] + (int128_t)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes: pa^2 mod(E)
//~ Note: E(X) = X^9 - 2
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128_t tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128_t)pa[0] * pa[0] + (((int128_t)pa[5] * pa[4] + (int128_t)pa[6] * pa[3] + (int128_t)pa[7] * pa[2] + (int128_t)pa[8] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128_t)pa[1] * pa[0]) << 1) + (((((int128_t)pa[6] * pa[4] + (int128_t)pa[7] * pa[3] + (int128_t)pa[8] * pa[2]) << 1) + (int128_t)pa[5] * pa[5]) << 1);
	tmp_prod_result[2] = (((int128_t)pa[2] * pa[0]) << 1) + (int128_t)pa[1] * pa[1] + (((int128_t)pa[6] * pa[5] + (int128_t)pa[7] * pa[4] + (int128_t)pa[8] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128_t)pa[2] * pa[1] + (int128_t)pa[3] * pa[0]) << 1) + (((((int128_t)pa[7] * pa[5] + (int128_t)pa[8] * pa[4]) << 1) + (int128_t)pa[6] * pa[6]) << 1);
	tmp_prod_result[4] = (((int128_t)pa[3] * pa[1] + (int128_t)pa[4] * pa[0]) << 1) + (int128_t)pa[2] * pa[2] + (((int128_t)pa[7] * pa[6] + (int128_t)pa[8] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128_t)pa[3] * pa[2] + (int128_t)pa[4] * pa[1] + (int128_t)pa[5] * pa[0]) << 1) + (((((int128_t)pa[8] * pa[6]) << 1) + (int128_t)pa[7] * pa[7]) << 1);
	tmp_prod_result[6] = (((int128_t)pa[4] * pa[2] + (int128_t)pa[5] * pa[1] + (int128_t)pa[6] * pa[0]) << 1) + (int128_t)pa[3] * pa[3] + (((int128_t)pa[8] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128_t)pa[4] * pa[3] + (int128_t)pa[5] * pa[2] + (int128_t)pa[6] * pa[1] + (int128_t)pa[7] * pa[0]) << 1) + (((int128_t)pa[8] * pa[8]) << 1);
	tmp_prod_result[8] = (((int128_t)pa[5] * pa[3] + (int128_t)pa[6] * pa[2] + (int128_t)pa[7] * pa[1] + (int128_t)pa[8] * pa[0]) << 1) + (int128_t)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ --------------------------------------------------------------------

//~ performs the internal reduction on 'op' and PUTS the result in 'rop'
void internal_reduction(int64_t *rop, int128_t *op){

	int64_t tmpQ[9];
	int128_t tmpZero[9];

	//~ computation of : op*neginv_red_int_mat mod(mont_phi), with result coeffs centered
	tmpQ[0] = (int64_t)op[0] * 0x3b390f5b34c60d42L + (int64_t)op[1] * -0x734db9a337d56229L + (int64_t)op[2] * -0x2c99c6bd6e12941dL + (int64_t)op[3] * 0x50de6e425b00807eL + (int64_t)op[4] * -0xa05f8a4b2918950L + (int64_t)op[5] * -0x632e78876ba4131L + (int64_t)op[6] * -0x2b4c33b0db7d297fL + (int64_t)op[7] * 0x57c207ff18dec928L + (int64_t)op[8] * 0x36037c87f4cd849L;
	tmpQ[1] = (int64_t)op[0] * 0xf435678c0678106L + (int64_t)op[1] * 0x7414f8f0f4ad442dL + (int64_t)op[2] * -0x2a88658e80cf2aa2L + (int64_t)op[3] * 0x1ac2d97aec89d568L + (int64_t)op[4] * -0x3d3249eaa7b8c35aL + (int64_t)op[5] * -0x7912d84c385cbba4L + (int64_t)op[6] * 0x7dd99accf0834bacL + (int64_t)op[7] * -0xdd7c1b9afcf048bL + (int64_t)op[8] * 0x28e5c8cfa54e5b5eL;
	tmpQ[2] = (int64_t)op[0] * 0x3428a69254e6009L + (int64_t)op[1] * -0x3c0c85b2779f58beL + (int64_t)op[2] * -0x4ab5c8c60621674eL + (int64_t)op[3] * -0x47485aadb69341dL + (int64_t)op[4] * -0x7763888e77483ceL + (int64_t)op[5] * 0x796e492d01c74c3L + (int64_t)op[6] * 0x42efa42c506366L + (int64_t)op[7] * -0x58efb720d4e5a6f9L + (int64_t)op[8] * 0x634a47f24efe6e45L;
	tmpQ[3] = (int64_t)op[0] * 0xe6fe0217d8e358aL + (int64_t)op[1] * -0x49ee89784949854aL + (int64_t)op[2] * -0x5b96b1f6ca4af2dcL + (int64_t)op[3] * 0x642b0f346015d218L + (int64_t)op[4] * -0x3fa15dfc8a2dd7e7L + (int64_t)op[5] * -0x7f02d03082c6996fL + (int64_t)op[6] * -0x6625004bfdf84ccL + (int64_t)op[7] * -0x52cb71c847f1cd1eL + (int64_t)op[8] * 0x32cb156a6f37f3b9L;
	tmpQ[4] = (int64_t)op[0] * 0xdaf762c45640468L + (int64_t)op[1] * -0x32fc16c22df218e9L + (int64_t)op[2] * -0x1ed0385ecf14ca0L + (int64_t)op[3] * 0x28c51cc5bca66a69L + (int64_t)op[4] * 0x16acd8ba846229e3L + (int64_t)op[5] * 0x40d6b95e91efcbedL + (int64_t)op[6] * -0x5f6b5bee7bb71c90L + (int64_t)op[7] * -0x3449a9cec0878387L + (int64_t)op[8] * -0x8d09d42aacf954aL;
	tmpQ[5] = (int64_t)op[0] * 0x3837dbbac6332f0aL + (int64_t)op[1] * -0x4dbef03d1a7bee51L + (int64_t)op[2] * 0x3b454664bac776baL + (int64_t)op[3] * -0x5e280aee0afcd9f3L + (int64_t)op[4] * -0x672cc2126c2121c9L + (int64_t)op[5] * 0x4eae7b1841bed078L + (int64_t)op[6] * 0x77aedb41defa8812L + (int64_t)op[7] * -0x52d056c041568593L + (int64_t)op[8] * 0xf0570db6bd5e25cL;
	tmpQ[6] = (int64_t)op[0] * -0x53bc4f898411621L + (int64_t)op[1] * 0x4fabf3c3076d3af1L + (int64_t)op[2] * -0x640a632a776de123L + (int64_t)op[3] * 0x477b3233869ab010L + (int64_t)op[4] * 0x265608b3da3155dcL + (int64_t)op[5] * 0x10bce0d639f84c18L + (int64_t)op[6] * -0x436531731e73048bL + (int64_t)op[7] * 0x5ba0f402ec261addL + (int64_t)op[8] * -0x2a645d33f69deb2cL;
	tmpQ[7] = (int64_t)op[0] * 0x1472e467d2a72dafL + (int64_t)op[1] * 0xf435678c0678106L + (int64_t)op[2] * 0x7414f8f0f4ad442dL + (int64_t)op[3] * -0x2a88658e80cf2aa2L + (int64_t)op[4] * 0x1ac2d97aec89d568L + (int64_t)op[5] * -0x3d3249eaa7b8c35aL + (int64_t)op[6] * -0x7912d84c385cbba4L + (int64_t)op[7] * 0x7dd99accf0834bacL + (int64_t)op[8] * -0xdd7c1b9afcf048bL;
	tmpQ[8] = (int64_t)op[0] * -0x787d47924a150ed2L + (int64_t)op[1] * 0x3837dbbac6332f0aL + (int64_t)op[2] * -0x4dbef03d1a7bee51L + (int64_t)op[3] * 0x3b454664bac776baL + (int64_t)op[4] * -0x5e280aee0afcd9f3L + (int64_t)op[5] * -0x672cc2126c2121c9L + (int64_t)op[6] * 0x4eae7b1841bed078L + (int64_t)op[7] * 0x77aedb41defa8812L + (int64_t)op[8] * -0x52d056c041568593L;

	//~ computation of : tmp_q*red_int_mat
	tmpZero[0] = (int128_t)tmpQ[0] * 0xb11ce988215efL + (int128_t)tmpQ[1] * -0x18505c68e1dfffL + (int128_t)tmpQ[2] * -0xe1d35454eb6e33L + (int128_t)tmpQ[3] * -0x150b35e6d810666L + (int128_t)tmpQ[4] * -0xf1f51e742e7ffL + (int128_t)tmpQ[5] * -0x11e7573f1f10572L + (int128_t)tmpQ[6] * -0x7881b92d7d05fdL + (int128_t)tmpQ[7] * -0xc947d8c4f5cb83L + (int128_t)tmpQ[8] * 0x2d294a9d41b3c9L;
	tmpZero[1] = (int128_t)tmpQ[0] * -0x1c282953d04abbL + (int128_t)tmpQ[1] * -0x5870c5d854e683L + (int128_t)tmpQ[2] * -0x19a04d25e9304eaL + (int128_t)tmpQ[3] * 0xf1f51e742e7ffL + (int128_t)tmpQ[4] * 0x2f30453c6b0779L + (int128_t)tmpQ[5] * 0xe1d35454eb6e33L + (int128_t)tmpQ[6] * -0xbc3f01586868c7L + (int128_t)tmpQ[7] * -0xa3588fd547c2edL + (int128_t)tmpQ[8] * -0xef452eb585fdf9L;
	tmpZero[2] = (int128_t)tmpQ[0] * -0x18fccdcbe69b18aL + (int128_t)tmpQ[1] * -0x715768cc78ac18L + (int128_t)tmpQ[2] * -0x91685c21a0cb92L + (int128_t)tmpQ[3] * -0x2f30453c6b0779L + (int128_t)tmpQ[4] * 0x1bb7534f84b536L + (int128_t)tmpQ[5] * 0x19a04d25e9304eaL + (int128_t)tmpQ[6] * 0x11ba492258cf616L + (int128_t)tmpQ[7] * -0x3f197e50e8ab24L + (int128_t)tmpQ[8] * 0xfd8aa7a4702369L;
	tmpZero[3] = (int128_t)tmpQ[0] * 0xf602f7965db2dbL + (int128_t)tmpQ[1] * 0x142fa3417244a8L + (int128_t)tmpQ[2] * -0x36b7f54a7d7807L + (int128_t)tmpQ[3] * -0x1bb7534f84b536L + (int128_t)tmpQ[4] * -0x1232c42339f66a4L + (int128_t)tmpQ[5] * 0x91685c21a0cb92L + (int128_t)tmpQ[6] * -0x119652297c1a8c9L + (int128_t)tmpQ[7] * 0xe6e3a5af1acadcL + (int128_t)tmpQ[8] * 0x76d8902af39e46L;
	tmpZero[4] = (int128_t)tmpQ[0] * 0x265c4e245bf9bbL + (int128_t)tmpQ[1] * -0x173a8843a370b2fL + (int128_t)tmpQ[2] * -0x6c3395561c59c1L + (int128_t)tmpQ[3] * 0x1232c42339f66a4L + (int128_t)tmpQ[4] * 0xc8e8e0912a043eL + (int128_t)tmpQ[5] * 0x36b7f54a7d7807L + (int128_t)tmpQ[6] * -0x62a4eaa8e09583L + (int128_t)tmpQ[7] * 0x558c9360c70134L + (int128_t)tmpQ[8] * 0x15a513cb2cacfd0L;
	tmpZero[5] = (int128_t)tmpQ[0] * -0x2ef2acd42f48eL + (int128_t)tmpQ[1] * -0x945786eee3c020L + (int128_t)tmpQ[2] * -0x3a14254c6e5dbaL + (int128_t)tmpQ[3] * -0xc8e8e0912a043eL + (int128_t)tmpQ[4] * -0xfa9b7cd72d2edbL + (int128_t)tmpQ[5] * 0x6c3395561c59c1L + (int128_t)tmpQ[6] * 0x121649ebf33fa00L + (int128_t)tmpQ[7] * 0x18c8288241c0a8L + (int128_t)tmpQ[8] * -0xc3e3878cafb6d4L;
	tmpZero[6] = (int128_t)tmpQ[0] * -0xcabbc7c37ad8dL + (int128_t)tmpQ[1] * -0x4363b1c6b52594L + (int128_t)tmpQ[2] * -0x11b11580b1faf7L + (int128_t)tmpQ[3] * 0xfa9b7cd72d2edbL + (int128_t)tmpQ[4] * -0x6a2cd3aebbf0b9L + (int128_t)tmpQ[5] * 0x3a14254c6e5dbaL + (int128_t)tmpQ[6] * -0xd52797d8b3c0a6L + (int128_t)tmpQ[7] * -0x12fd7feafd71431L + (int128_t)tmpQ[8] * 0x206c1a7606908L;
	tmpZero[7] = (int128_t)tmpQ[0] * -0x116b909720f13b7L + (int128_t)tmpQ[1] * -0x182ec9ec82b6d78L + (int128_t)tmpQ[2] * -0x1e244e42424de4L + (int128_t)tmpQ[3] * 0x6a2cd3aebbf0b9L + (int128_t)tmpQ[4] * -0x43983340bbbc19L + (int128_t)tmpQ[5] * 0x11b11580b1faf7L + (int128_t)tmpQ[6] * -0x834bc8ec83f133L + (int128_t)tmpQ[7] * -0x4dd028e0e50f79L + (int128_t)tmpQ[8] * -0x9840df44d5e5fL;
	tmpZero[8] = (int128_t)tmpQ[0] * 0x1b89d855ef33a1L + (int128_t)tmpQ[1] * -0x1e8a4cf67f2a19L + (int128_t)tmpQ[2] * 0x8f3ab9f8f882b9L + (int128_t)tmpQ[3] * 0x43983340bbbc19L + (int128_t)tmpQ[4] * 0xa859af36c08333L + (int128_t)tmpQ[5] * 0x1e244e42424de4L + (int128_t)tmpQ[6] * -0xacf234778fff33L + (int128_t)tmpQ[7] * -0xdf11a4813dfb3aL + (int128_t)tmpQ[8] * 0xba0ac4b7727e2aL;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
	rop[5] = (op[5] + tmpZero[5]) >> PHI_LOG2;
	rop[6] = (op[6] + tmpZero[6]) >> PHI_LOG2;
	rop[7] = (op[7] + tmpZero[7]) >> PHI_LOG2;
	rop[8] = (op[8] + tmpZero[8]) >> PHI_LOG2;
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

