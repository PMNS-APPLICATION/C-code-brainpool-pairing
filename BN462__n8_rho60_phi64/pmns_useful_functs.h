#ifndef USEFUL_FUNCTS
#define USEFUL_FUNCTS


void init_data();

void free_data();

void from_int_to_pmns(int64_t *rop, mpz_t op);

void from_pmns_to_int(mpz_t rop, int64_t *op);

void copy_poly(int64_t *rop, int64_t *op);

void print_element(int64_t *poly);

int cmp_poly_evals(int64_t *pa, int64_t *pb);

#endif

