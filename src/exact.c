#include "common.h"
#include <stdlib.h>

/* Exact Jaccard index, used as the "ground truth" to benchmark the sketch.
 * It keeps every canonical k-mer of both genomes in memory, so it is slow and
 * uses a lot of RAM on big genomes (about 8 bytes per k-mer). */

typedef struct {
    uint64_t *v;
    size_t n, cap;
    int failed;
} KmerBuf;

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void push_cb(uint64_t kmer, void *ctx) {
    KmerBuf *b = ctx;
    if (b->failed) return;
    if (b->n == b->cap) {
        size_t nc = b->cap ? b->cap * 2 : 1u << 16;
        uint64_t *p = realloc(b->v, nc * sizeof *p);
        if (!p) { b->failed = 1; return; }
        b->v = p;
        b->cap = nc;
    }
    b->v[b->n++] = kmer;
}

/* Collect, sort and de-duplicate. Returns the number of distinct k-mers,
 * or (size_t)-1 on out-of-memory. */
static size_t distinct_kmers(const Genome *g, int k, KmerBuf *b) {
    for_each_kmer(g, k, push_cb, b);
    if (b->failed) return (size_t)-1;
    qsort(b->v, b->n, sizeof *b->v, cmp_u64);
    size_t u = 0;
    for (size_t i = 0; i < b->n; i++)
        if (u == 0 || b->v[i] != b->v[u - 1]) b->v[u++] = b->v[i];
    b->n = u;
    return u;
}

/* J = shared / (|A| + |B| - shared). Returns 0.0 on error or if both are empty. */
double exact_jaccard(const Genome *a, const Genome *b, int k) {
    if (!a || !b || k < 1 || k > 31) return 0.0;

    KmerBuf A = {0}, B = {0};
    size_t na = distinct_kmers(a, k, &A);
    size_t nb = distinct_kmers(b, k, &B);
    double result = 0.0;

    if (na != (size_t)-1 && nb != (size_t)-1) {
        size_t i = 0, j = 0, shared = 0;
        while (i < na && j < nb) {
            if (A.v[i] == B.v[j])      { shared++; i++; j++; }
            else if (A.v[i] < B.v[j])  i++;
            else                       j++;
        }
        size_t uni = na + nb - shared;
        if (uni > 0) result = (double)shared / (double)uni;
    }
    free(A.v);
    free(B.v);
    return result;
}