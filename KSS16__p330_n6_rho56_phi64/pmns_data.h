#ifndef STRUCTS_DATA
#define STRUCTS_DATA


typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 5
#define NB_COEFF 6
#define NB_ADD_MAX 2

#define RHO_UPLOG2 56	// note: so 57 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 55

#define CONV_MASK 36028797018963967UL  // = (1 << CONVBASE_LOG2) - 1, for conversion



//Representations of polynomials Pi, used for conversion into the PMNS
//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x1e26e94d35eb6L, -0x20f0ca3cdd3122L, 0x34e127cb484d15L, -0x300a212956bbfdL, -0x9843ef706ad56L, -0x120c320ac3a295L},
	{0x29be3ec8ec3364L, -0x1df659bdbec165L, -0x7e559e38acb9eL, 0x58f0971ef7175L, 0x6806ca78feff7L, 0x2e63faa0828f27L},
	{-0x3e65c5cf2879bL, 0x1e87af6b07c480L, 0x13c26c19f5eb70L, -0x1c8a904957e891L, -0x20a671464b5861L, 0x4d4fb7e91d7b0eL},
	{-0x130fa9120630b0L, 0x2638feb3ce5295L, -0x48dd29585067d2L, -0x2d71c295e9fe51L, 0x1a79320843ee08L, 0x2dc6610a27b557L},
	{-0x116bbaddc1210aL, -0x1ab8c4cc80e62L, -0x914bf655fdf99L, 0x197d39b126ade7L, -0x49fd16f4165cf6L, -0xf91c522e35d83L},
	{0xde2d101519ab5L, 0x144203b414415L, -0x7a338cc8651f6L, 0x431bb0962ccd69L, 0x19ae9248606936L, -0x35e0c2ad96db5L}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

