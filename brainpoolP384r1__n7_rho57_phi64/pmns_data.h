#ifndef STRUCTS_DATA
#define STRUCTS_DATA


typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 6
#define NB_COEFF 7
#define NB_ADD_MAX 0

#define RHO_UPLOG2 57	// note: so 58 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 55

#define CONV_MASK 36028797018963967UL  // = (1 << CONVBASE_LOG2) - 1, for conversion



//Representations of polynomials Pi, used for conversion into the PMNS
//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{-0xa89e2e68de1c4L, 0xcc336f7c8c08dL, -0x2d9f83f53295bbL, -0x183cfcc550ce9bL, 0x9b7217e70b004L, 0xdfa98727f6ab2L, -0xd6c784265cb43L},
	{-0x427052b060402L, -0x629679161413f0L, -0x7d7699c9cf6d8L, 0x18b303c0d1f78eL, 0xfa72577b43ccfL, -0x25d87b6fbb2a93L, 0x26ab1f9e66c72aL},
	{-0x2b1e80107236cdL, 0x13c0947eb1f633L, 0x3d8698f1ecac57L, -0x2fe2a2f09dc39L, -0x1d118d60501bddL, 0x4044cbb68404c6L, 0x37786ae11e9b8bL},
	{-0x1f5e5d044588ccL, 0x1f40aba9acff84L, -0x1cccc0d843e244L, -0x32bfd0c029b10L, -0x92710c7ed4da7L, -0x2eeab32056f8bfL, 0x29c57f58630ea3L},
	{0x1ccabbc20be253L, -0x1fcf4e218634f6L, 0x2e0b6cbe0e0165L, 0x1214a98b34399fL, -0x200905f64a7433L, 0x1810374d71e569L, 0x180c12ecb58470L},
	{0xb7ef97b35f41cL, -0x3114950474383L, -0xc47e3b44fdd66L, 0x1dce83f2f782fcL, -0x16e45d3c27103cL, -0x1e3ac9b972f392L, -0x763f4645f2207L},
	{-0xfb3256ef78c88L, -0x127a0915c10425L, -0x1e781ffed67779L, -0x9c8ac3da659e5L, -0x108d15a892dcf5L, 0x14dbd35c1d2043L, -0x2b70cd804abe8fL}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

