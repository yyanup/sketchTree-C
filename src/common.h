#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>
#include <stddef.h>

/* ---- P1 provides ---- */
typedef struct { char name[64]; char *seq; size_t len; } Genome;
Genome *read_fasta(const char *path);              /* joins all records into one Genome */
void    free_genome(Genome *g);
/* calls cb(canonical_kmer) for every valid canonical k-mer */
void    for_each_kmer(const Genome *g, int k, void (*cb)(uint64_t, void *), void *ctx);

/* ---- P2 provides ---- */
typedef struct { char name[64]; int k; int n; uint64_t *h; } Sketch;  /* h sorted, n <= s */
uint64_t hash64(uint64_t x);
Sketch  *build_sketch(const Genome *g, int k, int s);
void     save_sketch(const Sketch *sk, const char *path);
Sketch  *load_sketch(const char *path);
double   sketch_jaccard(const Sketch *a, const Sketch *b, int s);
double   mash_distance(double j, int k);
double   exact_jaccard(const Genome *a, const Genome *b, int k);

/* ---- P3 provides ---- */
void upgma_newick(double **D, char **names, int n, const char *outfile);
#endif
