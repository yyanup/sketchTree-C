/* upgma.c - UPGMA tree building, Newick output (P3) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Newick-unsafe characters in a name become '_' */
static char *clean_name(const char *s) {
    char *r = malloc(strlen(s) + 1);
    if (!r) return NULL;
    size_t i = 0;
    for (; s[i]; i++) {
        char c = s[i];
        r[i] = (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' ||
                c == '[' || c == ']' || c == '\'' || c == ' ' || c == '\t' ||
                c == '\n' || c == '\r') ? '_' : c;
    }
    r[i] = '\0';
    return r;
}

void upgma_newick(double **D, char **names, int n, const char *outfile) {
    if (n < 1 || !D || !names || !outfile) {
        fprintf(stderr, "upgma_newick: bad arguments\n");
        return;
    }
    char **tree = calloc(n, sizeof(char *)); /* Newick text of each cluster */
    double *height = calloc(n, sizeof(double));
    int *size = malloc(n * sizeof(int));
    int *active = malloc(n * sizeof(int));
    double *d = malloc((size_t)n * n * sizeof(double)); /* private copy of D */
    if (!tree || !height || !size || !active || !d) {
        fprintf(stderr, "upgma_newick: out of memory\n");
        goto cleanup;
    }
    for (int i = 0; i < n; i++) {
        tree[i] = clean_name(names[i]);
        if (!tree[i]) { fprintf(stderr, "upgma_newick: out of memory\n"); goto cleanup; }
        size[i] = 1;
        active[i] = 1;
        for (int j = 0; j < n; j++) d[(size_t)i * n + j] = D[i][j];
    }

    for (int step = 0; step < n - 1; step++) {
        /* closest pair of active clusters */
        int a = -1, b = -1;
        double best = 0;
        for (int i = 0; i < n; i++) {
            if (!active[i]) continue;
            for (int j = i + 1; j < n; j++) {
                if (!active[j]) continue;
                double v = d[(size_t)i * n + j];
                if (a < 0 || v < best) { best = v; a = i; b = j; }
            }
        }
        double h = best / 2.0;
        double la = h - height[a], lb = h - height[b];
        if (la < 0) la = 0;   /* guard against non-ultrametric noise */
        if (lb < 0) lb = 0;

        size_t len = strlen(tree[a]) + strlen(tree[b]) + 80;
        char *node = malloc(len);
        if (!node) { fprintf(stderr, "upgma_newick: out of memory\n"); goto cleanup; }
        snprintf(node, len, "(%s:%.6f,%s:%.6f)", tree[a], la, tree[b], lb);
        free(tree[a]); free(tree[b]);
        tree[a] = node; tree[b] = NULL;

        /* distance from merged cluster to every other cluster */
        for (int x = 0; x < n; x++) {
            if (!active[x] || x == a || x == b) continue;
            double nd = (d[(size_t)a * n + x] * size[a] +
                         d[(size_t)b * n + x] * size[b]) / (double)(size[a] + size[b]);
            d[(size_t)a * n + x] = d[(size_t)x * n + a] = nd;
        }
        size[a] += size[b];
        height[a] = h;
        active[b] = 0;
    }

    int root = 0;
    while (root < n && !active[root]) root++;
    FILE *f = fopen(outfile, "w");
    if (!f) { fprintf(stderr, "upgma_newick: cannot write %s\n", outfile); goto cleanup; }
    fprintf(f, "%s;\n", tree[root]);
    fclose(f);

cleanup:
    if (tree) { for (int i = 0; i < n; i++) free(tree[i]); }
    free(tree); free(height); free(size); free(active); free(d);
}