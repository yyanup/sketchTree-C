#include "common.h"

/* A=0, C=1, G=2, T=3, anything else = -1 (ambiguous) */
static int base_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

/* Calls cb(canonical_kmer, ctx) for every valid k-mer. k must be 1..31. */
void for_each_kmer(const Genome *g, int k, void (*cb)(uint64_t, void *), void *ctx) {
    if (!g || !g->seq || !cb || k < 1 || k > 31 || g->len < (size_t)k) return;

    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;                            /* good bases in a row */

    for (size_t i = 0; i < g->len; i++) {
        int c = base_code(g->seq[i]);
        if (c < 0) {                          /* ambiguous base: start over */
            valid = 0;
            fwd = rev = 0;
            continue;
        }
        fwd = ((fwd << 2) | (uint64_t)c) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
        if (valid < k) valid++;
        if (valid >= k)
            cb(fwd < rev ? fwd : rev, ctx);   /* canonical = smaller of the two */
    }
}