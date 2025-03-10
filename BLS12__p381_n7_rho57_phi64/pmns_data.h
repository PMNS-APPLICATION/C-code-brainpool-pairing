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
	{-0x3b0d1a9fbcf6eL, 0x2445635b0b5171L, -0xe6bf0a3c7faaL, 0x3d971676b2718cL, -0x915aef5a8bb8cL, -0x79097a44bc48cL, 0x1a6b89371699a4L},
	{0x12112057a0288bL, -0x17b706a1ccce4cL, 0x8864577b2fa62L, 0x1bfcc97311c31aL, -0x33aa53cb470741L, 0x137bd9be6c464dL, 0x32947a6902cceL},
	{0x74f9fc1dde394L, -0x18f24d88d68749L, -0x407742985fa8fL, -0x31100acd5670eaL, -0x2bf91518744201L, 0xacd6984157e9bL, 0x33951733cf5e2L},
	{0x89d48f00b1490L, -0x2522a0bb2289abL, -0xbe8de5aa7e83fL, 0x8f903aaf85b06L, -0x38af3f399330fcL, 0x1934462ec72041L, -0x9e966b0aaba36L},
	{0x2bc7d096778223L, -0x1aac02a1bf1d08L, -0x28e9ede2bc03a4L, -0x1d10bcd38cf1d9L, -0x2f5b8b04b6aebeL, -0xd39d71db3802L, -0x2898106a3967cL},
	{0x2feb87fcdc6ec6L, -0x4de092fe687a9L, -0x328ac980fd8b2cL, 0xbc197cbbdac8bL, 0xec8b1317441f9L, -0x39332a81f3657L, 0x44a15a48c6890L},
	{-0x9e64ea8b0e668L, -0x160ed22a3fdafL, -0x1bf200cdf83936L, -0x39a0035048becdL, 0xa61e3da9042a3L, -0x3f4bb3e965aeL, -0x1d417237e335c4L}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

