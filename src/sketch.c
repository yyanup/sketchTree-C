#include "common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------- growable array of hashes ---------- */

typedef struct {
    uint64_t *v;
    size_t n, cap;
    int failed;                 /* set if realloc ever fails */
} HashBuf;

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void collect_cb(uint64_t kmer, void *ctx) {
    HashBuf *b = ctx;
    if (b->failed) return;
    if (b->n == b->cap) {
        size_t nc = b->cap ? b->cap * 2 : 1u << 16;
        uint64_t *p = realloc(b->v, nc * sizeof *p);
        if (!p) { b->failed = 1; return; }
        b->v = p;
        b->cap = nc;
    }
    b->v[b->n++] = hash64(kmer);
}

/* ---------- build_sketch ---------- */

/* Bottom-s MinHash: hash every canonical k-mer, sort, drop duplicates, keep
 * the s smallest distinct values. Returns NULL on failure. */
Sketch *build_sketch(const Genome *g, int k, int s) {
    if (!g || k < 1 || k > 31 || s < 1) return NULL;

    HashBuf b = {0};
    for_each_kmer(g, k, collect_cb, &b);
    if (b.failed || b.n == 0) {
        free(b.v);
        return NULL;                    /* out of memory, or genome shorter than k */
    }

    qsort(b.v, b.n, sizeof *b.v, cmp_u64);

    size_t u = 0;                       /* unique count, done in place */
    for (size_t i = 0; i < b.n; i++)
        if (u == 0 || b.v[i] != b.v[u - 1]) b.v[u++] = b.v[i];

    size_t keep = u < (size_t)s ? u : (size_t)s;

    Sketch *sk = malloc(sizeof *sk);
    uint64_t *h = malloc(keep * sizeof *h);   /* keep >= 1 here */
    if (!sk || !h) {
        free(sk);
        free(h);
        free(b.v);
        return NULL;
    }
    memcpy(h, b.v, keep * sizeof *h);
    free(b.v);

    strncpy(sk->name, g->name, sizeof sk->name - 1);
    sk->name[sizeof sk->name - 1] = '\0';
    sk->k = k;
    sk->n = (int)keep;
    sk->h = h;
    return sk;
}

/* ---------- save / load ---------- */

/* File layout: name[64], k (int), n (int), then n uint64 hash values. */
void save_sketch(const Sketch *sk, const char *path) {
    if (!sk || !path) return;
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "Cannot write %s (does the folder exist?)\n", path);
        return;
    }
    int ok = fwrite(sk->name, sizeof sk->name, 1, f) == 1
          && fwrite(&sk->k, sizeof sk->k, 1, f) == 1
          && fwrite(&sk->n, sizeof sk->n, 1, f) == 1
          && (sk->n == 0 || fwrite(sk->h, sizeof *sk->h, (size_t)sk->n, f) == (size_t)sk->n);
    if (!ok) fprintf(stderr, "Error while writing %s\n", path);
    fclose(f);
}

/* Returns NULL if the file is missing or damaged. Caller frees h and the Sketch. */
Sketch *load_sketch(const char *path) {
    if (!path) return NULL;
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", path);
        return NULL;
    }
    Sketch *sk = malloc(sizeof *sk);
    if (!sk) { fclose(f); return NULL; }
    sk->h = NULL;

    if (fread(sk->name, sizeof sk->name, 1, f) != 1 ||
        fread(&sk->k, sizeof sk->k, 1, f) != 1 ||
        fread(&sk->n, sizeof sk->n, 1, f) != 1 ||
        sk->k < 1 || sk->k > 31 || sk->n < 0) {
        fprintf(stderr, "%s is not a valid sketch file\n", path);
        fclose(f);
        free(sk);
        return NULL;
    }
    sk->name[sizeof sk->name - 1] = '\0';

    sk->h = malloc((sk->n > 0 ? (size_t)sk->n : 1) * sizeof *sk->h);
    if (!sk->h || (sk->n > 0 &&
        fread(sk->h, sizeof *sk->h, (size_t)sk->n, f) != (size_t)sk->n)) {
        fprintf(stderr, "%s is truncated\n", path);
        fclose(f);
        free(sk->h);
        free(sk);
        return NULL;
    }
    fclose(f);
    return sk;
}

/* ---------- Jaccard estimate and Mash distance ---------- */

/* Walk both sorted lists together and look at the s smallest values of their
 * union. j = (how many of those appear in both sketches) / (how many we looked
 * at). Stopping after s values is what keeps the estimate unbiased when the two
 * sketches are different sizes. */
double sketch_jaccard(const Sketch *a, const Sketch *b, int s) {
    if (!a || !b || s < 1 || a->k != b->k) return 0.0;

    int i = 0, j = 0, seen = 0, shared = 0;
    while (seen < s && i < a->n && j < b->n) {
        if (a->h[i] == b->h[j])      { shared++; i++; j++; }
        else if (a->h[i] < b->h[j])  i++;
        else                         j++;
        seen++;
    }
    /* one list ran out: the rest of the other still counts as union members */
    if (seen < s) {
        int left = (i < a->n) ? a->n - i : b->n - j;
        if (left > s - seen) left = s - seen;
        seen += left;
    }
    return seen > 0 ? (double)shared / seen : 0.0;
}

/* D = -ln(2j / (1 + j)) / k, the estimated per-base mutation rate. */
double mash_distance(double j, int k) {
    if (j <= 0.0) return 1.0;           /* nothing shared: maximally distant */
    if (j >= 1.0) return 0.0;           /* identical sketches */
    double d = -log(2.0 * j / (1.0 + j)) / (double)k;
    return d > 1.0 ? 1.0 : d;
}