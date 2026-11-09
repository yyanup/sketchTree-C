#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Reads a FASTA file and joins all records into one Genome. */
Genome *read_fasta(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "read_fasta: cannot open %s\n", path);
        return NULL;
    }

    Genome *g = calloc(1, sizeof(Genome));
    size_t cap = 1 << 20, len = 0;           /* start with 1 MB, double when full */
    char *buf = malloc(cap);
    if (!g || !buf) {
        fprintf(stderr, "read_fasta: out of memory\n");
        free(g); free(buf); fclose(f);
        return NULL;
    }

    char line[65536];
    int at_line_start = 1;                    /* are we at the start of a line? */
    int in_header = 0;                        /* is the current line a '>' header? */
    int have_name = 0;

    while (fgets(line, sizeof line, f)) {
        size_t n = strlen(line);
        if (at_line_start) in_header = (line[0] == '>');
        at_line_start = (n > 0 && line[n - 1] == '\n');

        if (in_header) {
            if (!have_name) {                 /* name = first word after '>' */
                size_t j = 0;
                for (size_t i = 1; i < n && j < sizeof g->name - 1; i++) {
                    if (isspace((unsigned char)line[i])) break;
                    g->name[j++] = line[i];
                }
                g->name[j] = '\0';
                have_name = 1;
            }
            continue;                         /* skip header text */
        }

        for (size_t i = 0; i < n; i++) {
            unsigned char c = (unsigned char)line[i];
            if (isspace(c)) continue;         /* skip newlines and spaces */
            if (len + 1 >= cap) {             /* grow the buffer */
                cap *= 2;
                char *tmp = realloc(buf, cap);
                if (!tmp) {
                    fprintf(stderr, "read_fasta: out of memory\n");
                    free(buf); free(g); fclose(f);
                    return NULL;
                }
                buf = tmp;
            }
            buf[len++] = (char)toupper(c);
        }
    }
    fclose(f);
    buf[len] = '\0';

    if (!have_name) {                         /* no header: use the file name */
        const char *base = path, *p;
        for (p = path; *p; p++)
            if (*p == '/' || *p == '\\') base = p + 1;
        strncpy(g->name, base, sizeof g->name - 1);
    }
    g->seq = buf;
    g->len = len;
    return g;
}

void free_genome(Genome *g) {
    if (!g) return;
    free(g->seq);
    free(g);
}