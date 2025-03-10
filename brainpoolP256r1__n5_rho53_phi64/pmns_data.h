#ifndef STRUCTS_DATA
#define STRUCTS_DATA


typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 4
#define NB_COEFF 5
#define NB_ADD_MAX 1

#define RHO_UPLOG2 53	// note: so 54 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 52

#define CONV_MASK 4503599627370495UL  // = (1 << CONVBASE_LOG2) - 1, for conversion



//Representations of polynomials Pi, used for conversion into the PMNS
//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x10106a9110360L, 0xdbae7a042939L, 0x3fa5ba2640225L, 0x27a1866edb270L, -0x1482394641618L},
	{-0xcf667941bf13L, -0x533ff47ae6a49L, 0x3bc28b57d93c1L, 0x5276ae779b08L, 0x16d1a281e9869L},
	{0x3a5083bd8a544L, 0x2e2427ba4991aL, -0x2986b7f6a5b41L, -0x1460d22ac17fcL, 0x1e4e10f410a4fL},
	{-0x2739d2c00b24cL, 0x2be174f95d9a6L, -0x22f5f3c34ac8cL, 0x3e11370971c1L, 0x54b7b12494ffL},
	{0x2a35ddfc6c304L, -0x146dd1952e632L, -0x21eacb785c4c4L, -0x15f9d7c4c6c46L, -0x5ed626806d71eL}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

