/* Temporary test for P1's code. NOT part of src/, so it can have its own main().
   Build: gcc -std=c11 -Wall -Wextra -g -o test_p1 test_p1.c src/fasta.c src/kmer.c -lm
   Run:   .\test_p1.exe                      (built-in k-mer tests)
          .\test_p1.exe data\lambda.fna      (also tests the FASTA reader) */
#include <stdio.h>
#include <stdlib.h>
#include "src/common.h"

static void print_kmer(uint64_t x, void *ctx) {
    (void)ctx;
    printf("  canonical = %llu\n", (unsigned long long)x);
}
static void count_kmer(uint64_t x, void *ctx) {
    (void)x;
    (*(unsigned long long *)ctx)++;
}

int main(int argc, char **argv) {
    char s1[] = "ACGT";
    Genome t1 = {"t1", s1, 4};
    printf("Test 1: ACGT, k=3 (expect 2 k-mers, both canonical = 6)\n");
    for_each_kmer(&t1, 3, print_kmer, NULL);

    char s2[] = "ACNGT";
    Genome t2 = {"t2", s2, 5};
    unsigned long long n = 0;
    for_each_kmer(&t2, 2, count_kmer, &n);
    printf("Test 2: ACNGT, k=2 (expect 2 k-mers, skips the N): got %llu\n", n);

    if (argc > 1) {
        Genome *g = read_fasta(argv[1]);
        if (!g) return 1;
        n = 0;
        for_each_kmer(g, 21, count_kmer, &n);
        printf("Test 3: %s\n  name = %s\n  length = %zu bases (compare with NCBI)\n  k=21 k-mers = %llu\n",
               argv[1], g->name, g->len, n);
        free_genome(g);
    }
    return 0;
}