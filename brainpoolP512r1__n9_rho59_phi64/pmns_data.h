#ifndef STRUCTS_DATA
#define STRUCTS_DATA


typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 8
#define NB_COEFF 9
#define NB_ADD_MAX 0

#define RHO_UPLOG2 59	// note: so 60 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 57

#define CONV_MASK 144115188075855871UL  // = (1 << CONVBASE_LOG2) - 1, for conversion



//Representations of polynomials Pi, used for conversion into the PMNS
//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x1e4e9cf8539316L, 0xa6d828c52a1864L, 0x30f37cad404201L, 0x1c3208a1393abbL, -0xa9de2703d7ae59L, 0x5a7f3285ede1c7L, -0x143c83ee6c5120L, 0x2ccaf957ae5109L, 0x115e7400f5e790L},
	{-0x4992592780b8e8L, -0xd2899fd7938dffL, 0x1edce28c7fdf02L, 0xa7cddf1306cc26L, 0x6802ab7becf57fL, 0x885e16598a5ca7L, -0x957b687b5daee9L, -0x28bddf3597984dL, -0x992081468f6ecbL},
	{0x53bc54d11a7b37L, 0x33fff0d85326f6L, 0x694c7cb284c12aL, -0x7c3f1d371f2ad0L, 0xb51b34a074db81L, -0x333d7d868eefd3L, -0x6c1fe39b4ffbd8L, 0xa3c9b6928cbb0bL, 0x592eadcadb78cL},
	{-0x12358c66c731fcbL, -0x835c742eaa9d7dL, 0x41080021222ad6L, 0xc72e5ff4dd337bL, 0x185877ee6c9e51eL, -0x860b339c336998L, -0x33a2f1bec84fcbL, -0xbc27babea8d42L, 0x4d337dd5c2348fL},
	{0x14e6254566c987L, 0xe71890aa8965e9L, -0x6cc57da976ed79L, -0x40f85dbe07722cL, -0x8a21ebb576be75L, -0x20bf033e360b43L, -0x3326ca18e4789bL, -0x13798502c044c8dL, -0x196fd9cea567dbL},
	{0x161b7d77497caeL, -0x8dd752aa806748L, 0xd0afd537395f5fL, -0x8f717ce0b8712dL, -0x4e73c2a18bbb14L, -0xe5991846602e1L, -0xd9497dd53bd2e7L, -0x48aeb86bfd8bcbL, 0x3febc97ac93aaL},
	{0x9b64a1e32596d9L, -0xb3d14e4cf58405L, -0x1f0b5c0e3b75fdL, 0x305ce04fdcedd7L, 0xd76cb5f56c9efL, -0x8f1eda74399f8dL, 0x3dffcfa6b5b2b6L, -0x1f776541d6a447L, 0x8fe8fe92eaa5aaL},
	{-0x484d38ebef3b9L, -0x1333e3703a88d22L, 0x5f20fad1ff0ec6L, 0x177c782d2b8732L, 0xc90b34cb50e00eL, 0x71821716f460e8L, -0x26b69efe34f272L, 0xe39c95c4789053L, -0x5a096670bae654L},
	{0x3a648e61ad2d3fL, -0xe76057055ade65L, -0x199fc608bd9ad4aL, 0xe6a2219074e3e5L, 0xa7f7af0ba5c689L, -0xd1d0e47544b0aeL, -0x95eaa9d1f098d5L, -0x64fe6b0d4c9cb4L, 0x88912b6e8c69baL}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

