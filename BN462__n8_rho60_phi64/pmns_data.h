#ifndef STRUCTS_DATA
#define STRUCTS_DATA


typedef __int128 int128_t;

#define PHI_LOG2 64
#define POLY_DEG 7
#define NB_COEFF 8
#define NB_ADD_MAX 0

#define RHO_UPLOG2 60	// note: so 61 bits are used to store each coeficient of an element

#define CONVBASE_LOG2 58

#define CONV_MASK 288230376151711743UL  // = (1 << CONVBASE_LOG2) - 1, for conversion



//Representations of polynomials Pi, used for conversion into the PMNS
//Note that each Pi is a (random) representation of: (CONV_MASK+1)^i * phi^2
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{-0x1bdad347086273L, 0x579c04d1004742L, -0x118742062d6bb95L, -0x1299bf7b8e2c2f9L, 0x32b8688fd4db05L, 0xb6e28cf2ada767L, -0x442972fdc5507dL, -0x11373164686bea8L},
	{-0x6fac4ce2de0513L, 0x122e01778e64d50L, 0x903b4559f608a0L, -0x7e420eb48ce36dL, 0x483e0b660d0ae8L, -0x7f1683b08cb511L, -0x7bf0b59e6a91e6L, -0x201fb5a478fc20dL},
	{-0xbd37fe37d898adL, 0x139121604efc65cL, -0xe062531586155dL, -0x18ccff900677b8cL, 0x2a2d9488b1916bL, -0x3994205c3d368fL, 0x5ae657b847faa5L, -0x194d250f5bd9ab9L},
	{0x17fabf4fde9a6a9L, 0x145d7b7789a4139L, 0x17ca06dd559e05fL, -0xade81b3f48d81fL, -0x94dcfa423dbdafL, -0x132364c032c8c8fL, 0x9901c92c710623L, 0x1105bbc7a15f134L},
	{0x20710499c792782L, -0x5736c3917292fdL, 0xc75bb725858902L, 0x160fe8efae5eddL, -0x2d8b8fe1b1c824L, -0x31e11bb000f54eL, 0x1440a9261984fc8L, 0x119175e46bccbd4L},
	{-0x75cd42b5232a5eL, -0x70453e8954aa0dL, 0xbec1f0b0c6c648L, -0x997c9844840331L, 0x1a1f60ec7677dcbL, -0x10e54b95e1224d5L, 0xce91be031fd8c3L, -0xe4c027cbd98f64L},
	{-0x1aaba1a00525519L, 0x62e9df0ad6daf7L, -0x19d28346a4155e9L, 0x307c37c792958L, -0x21254358bb4841L, 0x267474774d23ef4L, 0xc2be07bc6f4bf3L, -0x4b6e60513ffb6dL},
	{-0x246154df7af8431L, -0x86926a411eb1ebL, 0x593e8d15285a23L, -0x4df00743bdb700L, -0x170f534f22fdbfcL, 0xd43b377e28a67L, -0x8f53e577920a99L, -0x1379611505d265fL}};


mpz_t modul_p;

mpz_t gama_pow[POLY_DEG];

#endif

