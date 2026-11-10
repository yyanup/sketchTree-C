#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_GENOMES 256
#define LINE_BUF    65536

static void usage(void) {
    fprintf(stderr,
        "SketchTree-C: fast genome distances and trees using MinHash sketches\n\n"
        "Usage:\n"
        "  sketchtree info   <fasta> <k>                      genome length and k-mer count\n"
        "  sketchtree sketch <fasta> <k> <s> <out.sk>         build and save one sketch\n"
        "  sketchtree dist   <taxonomy.csv> <data_dir> <k> <s> <out.csv>\n"
        "                                                     distance matrix of all genomes\n"
        "  sketchtree tree   <matrix.csv> <out.nwk>           UPGMA tree in Newick format\n\n"
        "k = k-mer length (1..31), s = sketch size (e.g. 1000)\n\n"
        "Example:\n"
        "  sketchtree dist data/taxonomy.csv data 21 1000 results/dist_sketch.csv\n"
        "  sketchtree tree results/dist_sketch.csv results/tree.nwk\n");
}

/* ---------- small helpers ---------- */

static int parse_int(const char *str, int lo, int hi, const char *what, int *out) {
    char *end;
    long v = strtol(str, &end, 10);
    if (*str == '\0' || *end != '\0' || v < lo || v > hi) {
        fprintf(stderr, "Invalid %s '%s' (must be a whole number from %d to %d)\n",
                what, str, lo, hi);
        return 0;
    }
    *out = (int)v;
    return 1;
}

static char *dup_str(const char *s) {
    size_t n = strlen(s) + 1;
    char *p = malloc(n);
    if (p) memcpy(p, s, n);
    return p;
}

static void trim_eol(char *s) {
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r')) s[--n] = '\0';
}

/* Splits a line at commas (in place). Returns the number of fields. */
static int split_csv(char *line, char **fields, int max) {
    int n = 0;
    char *p = line;
    while (n < max) {
        fields[n++] = p;
        char *c = strchr(p, ',');
        if (!c) break;
        *c = '\0';
        p = c + 1;
    }
    return n;
}

static void free_sketch(Sketch *sk) {
    if (!sk) return;
    free(sk->h);
    free(sk);
}

static void count_cb(uint64_t kmer, void *ctx) {
    (void)kmer;
    (*(unsigned long long *)ctx)++;
}

/* ---------- command: info ---------- */

static int cmd_info(const char *path, int k) {
    Genome *g = read_fasta(path);
    if (!g) return 1;
    unsigned long long n = 0;
    for_each_kmer(g, k, count_cb, &n);
    printf("File:    %s\nName:    %s\nLength:  %zu bases\nk-mers:  %llu (k=%d)\n",
           path, g->name, g->len, n, k);
    free_genome(g);
    return 0;
}

/* ---------- command: sketch ---------- */

static int cmd_sketch(const char *path, int k, int s, const char *out) {
    Genome *g = read_fasta(path);
    if (!g) return 1;
    Sketch *sk = build_sketch(g, k, s);
    free_genome(g);
    if (!sk) {
        fprintf(stderr, "Could not build a sketch for %s\n", path);
        return 1;
    }
    save_sketch(sk, out);
    printf("Sketch of %s saved to %s (%d hash values, k=%d)\n", path, out, sk->n, k);
    free_sketch(sk);
    return 0;
}

/* ---------- command: dist ---------- */

static int cmd_dist(const char *taxfile, const char *dir, int k, int s, const char *out) {
    FILE *f = fopen(taxfile, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", taxfile);
        return 1;
    }

    Sketch *sk[MAX_GENOMES];
    char   *labels[MAX_GENOMES];
    int n = 0, rc = 1;
    double **D = NULL;
    static char line[LINE_BUF];

    if (!fgets(line, sizeof line, f)) {       /* skip the header line */
        fprintf(stderr, "%s is empty\n", taxfile);
        fclose(f);
        return 1;
    }

    while (fgets(line, sizeof line, f)) {
        trim_eol(line);
        if (line[0] == '\0') continue;
        if (n >= MAX_GENOMES) {
            fprintf(stderr, "Too many genomes (max %d)\n", MAX_GENOMES);
            goto cleanup;
        }

        char *fields[8];
        split_csv(line, fields, 8);
        const char *file = fields[0];

        char label[64];
        strncpy(label, file, sizeof label - 1);
        label[sizeof label - 1] = '\0';
        char *dot = strrchr(label, '.');          /* ecoli.fna -> ecoli */
        if (dot) *dot = '\0';

        char path[512];
        snprintf(path, sizeof path, "%s/%s", dir, file);

        Genome *g = read_fasta(path);
        if (!g) goto cleanup;
        sk[n] = build_sketch(g, k, s);
        free_genome(g);
        if (!sk[n]) {
            fprintf(stderr, "Could not build a sketch for %s\n", path);
            goto cleanup;
        }
        labels[n] = dup_str(label);
        printf("[%d] %s: %d hash values\n", n + 1, label, sk[n]->n);
        n++;
    }

    if (n < 2) {
        fprintf(stderr, "Need at least 2 genomes, found %d\n", n);
        goto cleanup;
    }

    D = malloc(n * sizeof *D);
    if (!D) goto cleanup;
    for (int i = 0; i < n; i++) {
        D[i] = calloc(n, sizeof **D);
        if (!D[i]) goto cleanup;
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double jac = sketch_jaccard(sk[i], sk[j], s);
            double d = mash_distance(jac, k);
            D[i][j] = D[j][i] = d;
        }
    }

    FILE *o = fopen(out, "w");
    if (!o) {
        fprintf(stderr, "Cannot write %s (does the folder exist?)\n", out);
        goto cleanup;
    }
    fprintf(o, "name");
    for (int j = 0; j < n; j++) fprintf(o, ",%s", labels[j]);
    fprintf(o, "\n");
    for (int i = 0; i < n; i++) {
        fprintf(o, "%s", labels[i]);
        for (int j = 0; j < n; j++) fprintf(o, ",%.6f", D[i][j]);
        fprintf(o, "\n");
    }
    fclose(o);
    printf("Distance matrix (%d x %d) written to %s\n", n, n, out);
    rc = 0;

cleanup:
    fclose(f);
    if (D) {
        for (int i = 0; i < n; i++) free(D[i]);
        free(D);
    }
    for (int i = 0; i < n; i++) {
        free_sketch(sk[i]);
        free(labels[i]);
    }
    return rc;
}

/* ---------- command: tree ---------- */

static int cmd_tree(const char *csv, const char *out) {
    FILE *f = fopen(csv, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", csv);
        return 1;
    }

    static char line[LINE_BUF];
    char *fields[MAX_GENOMES + 2];
    char **names = NULL;
    double **D = NULL;
    int n = 0, rows = 0, rc = 1;

    if (!fgets(line, sizeof line, f)) {
        fprintf(stderr, "%s is empty\n", csv);
        goto cleanup;
    }
    trim_eol(line);
    int nf = split_csv(line, fields, MAX_GENOMES + 2);
    n = nf - 1;                                   /* first cell is just a corner label */
    if (n < 2) {
        fprintf(stderr, "%s needs at least 2 genomes in its header row\n", csv);
        n = 0;
        goto cleanup;
    }

    names = calloc(n, sizeof *names);
    D = calloc(n, sizeof *D);
    if (!names || !D) goto cleanup;
    for (int i = 0; i < n; i++) names[i] = dup_str(fields[i + 1]);

    while (rows < n && fgets(line, sizeof line, f)) {
        trim_eol(line);
        if (line[0] == '\0') continue;
        int c = split_csv(line, fields, MAX_GENOMES + 2);
        if (c != n + 1) {
            fprintf(stderr, "Row %d of %s has %d values, expected %d\n",
                    rows + 1, csv, c - 1, n);
            goto cleanup;
        }
        D[rows] = malloc(n * sizeof **D);
        if (!D[rows]) goto cleanup;
        for (int j = 0; j < n; j++) {
            char *end;
            D[rows][j] = strtod(fields[j + 1], &end);
            if (end == fields[j + 1]) {
                fprintf(stderr, "Bad number '%s' in row %d of %s\n",
                        fields[j + 1], rows + 1, csv);
                goto cleanup;
            }
        }
        rows++;
    }
    if (rows != n) {
        fprintf(stderr, "%s has %d data rows, expected %d\n", csv, rows, n);
        goto cleanup;
    }

    upgma_newick(D, names, n, out);
    printf("UPGMA tree of %d genomes written to %s\n", n, out);
    rc = 0;

cleanup:
    fclose(f);
    if (D) {
        for (int i = 0; i < n; i++) free(D[i]);
        free(D);
    }
    if (names) {
        for (int i = 0; i < n; i++) free(names[i]);
        free(names);
    }
    return rc;
}

/* ---------- main ---------- */

int main(int argc, char **argv) {
    if (argc < 2) {
        usage();
        return 1;
    }
    const char *cmd = argv[1];
    int k, s;

    if (strcmp(cmd, "info") == 0 && argc == 4) {
        if (!parse_int(argv[3], 1, 31, "k", &k)) return 1;
        return cmd_info(argv[2], k);
    }
    if (strcmp(cmd, "sketch") == 0 && argc == 6) {
        if (!parse_int(argv[3], 1, 31, "k", &k)) return 1;
        if (!parse_int(argv[4], 1, 10000000, "s", &s)) return 1;
        return cmd_sketch(argv[2], k, s, argv[5]);
    }
    if (strcmp(cmd, "dist") == 0 && argc == 7) {
        if (!parse_int(argv[4], 1, 31, "k", &k)) return 1;
        if (!parse_int(argv[5], 1, 10000000, "s", &s)) return 1;
        return cmd_dist(argv[2], argv[3], k, s, argv[6]);
    }
    if (strcmp(cmd, "tree") == 0 && argc == 4) {
        return cmd_tree(argv[2], argv[3]);
    }

    usage();
    return 1;
}