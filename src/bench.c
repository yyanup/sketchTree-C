/* bench.c - benchmark: sketch distance vs exact distance (P3)
 *
 *   int run_bench(const char *datadir, const char *taxfile, const char *outfile);
 *
 * Example: run_bench("data", "data/taxonomy.csv", "results/benchmark.csv");
 * Writes outfile (accuracy/time/memory per k,s) and a second file
 * "<outfile minus .csv>_time.csv" (per genome: length and sketching time).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "common.h"

#define MAXG 64
#define MAXF 32

/* ---------- small CSV helpers ---------- */
static char *trim(char *s) {
    while (*s == ' ' || *s == '\t' || *s == '"') s++;
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == ' ' || s[n-1] == '\t' || s[n-1] == '"' ||
                     s[n-1] == '\r' || s[n-1] == '\n')) s[--n] = '\0';
    return s;
}

static int split_csv(char *line, char **f, int maxf) {
    int n = 0;
    char *p = line;
    while (n < maxf) {
        f[n++] = p;
        char *c = strchr(p, ',');
        if (!c) break;
        *c = '\0';
        p = c + 1;
    }
    for (int i = 0; i < n; i++) f[i] = trim(f[i]);
    return n;
}

static int has_word(const char *hay, const char *needle) {
    char a[128]; size_t i;
    for (i = 0; hay[i] && i < sizeof a - 1; i++) a[i] = (char)tolower((unsigned char)hay[i]);
    a[i] = '\0';
    return strstr(a, needle) != NULL;
}

static int find_col(char **h, int n, const char *key) {
    for (int i = 0; i < n; i++) if (has_word(h[i], key)) return i;
    return -1;
}

static Genome *load_genome(const char *dir, const char *id) {
    static const char *ext[] = {"", ".fasta", ".fa", ".fna", ".fas"};
    char path[512];
    for (size_t i = 0; i < sizeof ext / sizeof ext[0]; i++) {
        snprintf(path, sizeof path, "%s/%s%s", dir, id, ext[i]);
        FILE *f = fopen(path, "rb");
        if (f) { fclose(f); return read_fasta(path); }
    }
    return NULL;
}

/* ---------- main benchmark ---------- */
int run_bench(const char *datadir, const char *taxfile, const char *outfile) {
    static Genome *G[MAXG];
    static char phylum[MAXG][64], genus[MAXG][64];
    char line[4096], *f[MAXF];
    int n = 0;

    FILE *tf = fopen(taxfile, "r");
    if (!tf) { fprintf(stderr, "bench: cannot open %s\n", taxfile); return 1; }
    if (!fgets(line, sizeof line, tf)) { fclose(tf); fprintf(stderr, "bench: empty taxonomy\n"); return 1; }
    int nf = split_csv(line, f, MAXF);

    int cfile = -1;
    const char *keys[] = {"file", "path", "accession", "assembly", "genome", "name"};
    for (size_t i = 0; i < sizeof keys / sizeof keys[0] && cfile < 0; i++)
        cfile = find_col(f, nf, keys[i]);
    if (cfile < 0) cfile = 0;
    int cphy = find_col(f, nf, "phylum");
    int cgen = find_col(f, nf, "genus");
    printf("bench: taxonomy columns -> file/name: '%s', phylum: %s, genus: %s\n",
           f[cfile], cphy >= 0 ? f[cphy] : "(none)", cgen >= 0 ? f[cgen] : "(none)");
    if (cphy < 0) fprintf(stderr, "bench: WARNING no 'phylum' column found; phylum accuracy will be 0\n");

    while (n < MAXG && fgets(line, sizeof line, tf)) {
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
        int c = split_csv(line, f, MAXF);
        if (cfile >= c) continue;
        Genome *g = load_genome(datadir, f[cfile]);
        if (!g) { fprintf(stderr, "bench: cannot find genome file for '%s' in %s (skipped)\n", f[cfile], datadir); continue; }
        G[n] = g;
        snprintf(phylum[n], 64, "%s", (cphy >= 0 && cphy < c) ? f[cphy] : "");
        snprintf(genus[n], 64, "%s", (cgen >= 0 && cgen < c) ? f[cgen] : "");
        n++;
    }
    fclose(tf);
    if (n < 2) { fprintf(stderr, "bench: need at least 2 genomes, loaded %d\n", n); return 1; }
    printf("bench: loaded %d genomes\n", n);

    FILE *out = fopen(outfile, "w");
    if (!out) { fprintf(stderr, "bench: cannot write %s\n", outfile); return 1; }
    fprintf(out, "k,s,pearson_r,mae,sketch_time_sec,mem_per_sketch_bytes,nn_phylum_acc,nn_genus_acc,pairs_used\n");

    char tpath[512];
    snprintf(tpath, sizeof tpath, "%s", outfile);
    size_t L = strlen(tpath);
    if (L > 4 && strcmp(tpath + L - 4, ".csv") == 0) tpath[L - 4] = '\0';
    strncat(tpath, "_time.csv", sizeof tpath - strlen(tpath) - 1);
    FILE *tout = fopen(tpath, "w");
    if (tout) fprintf(tout, "k,s,genome,genome_len,sketch_time_sec\n");

    const int ks[] = {15, 21, 27};
    const int ss[] = {100, 1000, 10000};
    double *ex = malloc((size_t)n * n * sizeof(double));
    double *sd = malloc((size_t)n * n * sizeof(double));
    if (!ex || !sd) { fprintf(stderr, "bench: out of memory\n"); return 1; }

    for (int ki = 0; ki < 3; ki++) {
        int k = ks[ki];
        /* exact distances: computed once per k */
        for (int i = 0; i < n; i++) {
            ex[(size_t)i * n + i] = 0;
            for (int j = i + 1; j < n; j++) {
                double d = mash_distance(exact_jaccard(G[i], G[j], k), k);
                ex[(size_t)i * n + j] = ex[(size_t)j * n + i] = d;
            }
        }
        printf("bench: k=%d exact distances done\n", k);

        for (int si = 0; si < 3; si++) {
            int s = ss[si];
            static Sketch *sk[MAXG];
            double total_t = 0, mem = 0;
            for (int i = 0; i < n; i++) {
                clock_t t0 = clock();
                sk[i] = build_sketch(G[i], k, s);
                double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
                total_t += dt;
                if (sk[i]) mem += (double)sk[i]->n * 8.0;
                if (tout) fprintf(tout, "%d,%d,%s,%zu,%.6f\n", k, s, G[i]->name, G[i]->len, dt);
            }
            mem /= n;

            for (int i = 0; i < n; i++) {
                sd[(size_t)i * n + i] = 0;
                for (int j = i + 1; j < n; j++) {
                    double d = mash_distance(sketch_jaccard(sk[i], sk[j], s), k);
                    sd[(size_t)i * n + j] = sd[(size_t)j * n + i] = d;
                }
            }

            /* Pearson r and MAE over all pairs (skip non-finite values) */
            double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, mae = 0;
            long m = 0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) {
                    double x = sd[(size_t)i * n + j], y = ex[(size_t)i * n + j];
                    if (!isfinite(x) || !isfinite(y)) continue;
                    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
                    mae += fabs(x - y);
                    m++;
                }
            double r = NAN;
            if (m > 1) {
                double cov = sxy - sx * sy / m;
                double vx = sxx - sx * sx / m, vy = syy - sy * sy / m;
                if (vx > 0 && vy > 0) r = cov / sqrt(vx * vy);
                mae /= m;
            } else mae = NAN;

            /* nearest neighbour check (using sketch distance) */
            int okp = 0, okg = 0;
            for (int i = 0; i < n; i++) {
                int best = -1;
                for (int j = 0; j < n; j++) {
                    if (j == i) continue;
                    double d = sd[(size_t)i * n + j];
                    if (!isfinite(d)) continue;
                    if (best < 0 || d < sd[(size_t)i * n + best]) best = j;
                }
                if (best < 0) continue;
                if (phylum[i][0] && strcmp(phylum[i], phylum[best]) == 0) okp++;
                if (genus[i][0] && strcmp(genus[i], genus[best]) == 0) okg++;
            }

            fprintf(out, "%d,%d,%.6f,%.6f,%.6f,%.0f,%.4f,%.4f,%ld\n",
                    k, s, r, mae, total_t, mem, (double)okp / n, (double)okg / n, m);
            printf("bench: k=%d s=%d  r=%.4f  MAE=%.5f  time=%.3fs  phylum-NN=%d/%d\n",
                   k, s, r, mae, total_t, okp, n);
            fflush(out);
            /* sketches are small; left allocated on purpose (no free_sketch in common.h) */
        }
    }
    fclose(out);
    if (tout) fclose(tout);
    free(ex); free(sd);
    for (int i = 0; i < n; i++) free_genome(G[i]);
    printf("bench: wrote %s and %s\n", outfile, tpath);
    return 0;
}