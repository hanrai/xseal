#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "minimap.h"
#include "kseq.h"
#include "mmpriv.h"

// KSEQ_INIT is needed for bseq.c
// But Minimap2 already has it.

typedef struct {
    int k, w;
    const char *fn;
    int n_threads;
    uint64_t total_bases;
    uint64_t total_hits;
} bench_args_t;

typedef struct {
    int k, w;
    int n_seqs;
    mm_bseq1_t *seqs;
    uint64_t hits;
} thread_data_t;

void *worker(void *data) {
    thread_data_t *d = (thread_data_t*)data;
    mm128_v v = {0,0,0};
    for (int i = 0; i < d->n_seqs; ++i) {
        v.n = 0;
        mm_sketch(NULL, d->seqs[i].seq, d->seqs[i].l_seq, d->w, d->k, i, 0, &v);
        d->hits += v.n;
    }
    free(v.a);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <fasta> [k] [w] [threads]\n", argv[0]);
        return 1;
    }

    int k = argc > 2 ? atoi(argv[2]) : 31;
    int w = argc > 3 ? atoi(argv[3]) : 21;
    int n_threads = argc > 4 ? atoi(argv[4]) : 12;
    const char *fn = argv[1];

    printf("🚀 Minimap2 Sketch Bench: mm_sketch\n");
    printf("--- Parameters: K=%d, W=%d, Threads=%d ---\n", k, w, n_threads);

    mm_bseq_file_t *fp = mm_bseq_open(fn);
    if (!fp) {
        perror("mm_bseq_open");
        return 1;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    uint64_t total_bases = 0;
    uint64_t total_hits = 0;

    // Process in batches to match Minimap2's style
    int batch_size = 100000000; // 100MB batch
    while (1) {
        int n_seqs;
        mm_bseq1_t *seqs = mm_bseq_read(fp, batch_size, 0, &n_seqs);
        if (!seqs) break;

        pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
        thread_data_t *tdata = malloc(n_threads * sizeof(thread_data_t));

        int seqs_per_thread = (n_seqs + n_threads - 1) / n_threads;
        for (int i = 0; i < n_threads; ++i) {
            tdata[i].k = k;
            tdata[i].w = w;
            tdata[i].hits = 0;
            tdata[i].n_seqs = (i == n_threads - 1) ? (n_seqs - i * seqs_per_thread) : seqs_per_thread;
            if (tdata[i].n_seqs < 0) tdata[i].n_seqs = 0;
            tdata[i].seqs = &seqs[i * seqs_per_thread];
            pthread_create(&threads[i], NULL, worker, &tdata[i]);
        }

        for (int i = 0; i < n_threads; ++i) {
            pthread_join(threads[i], NULL);
            total_hits += tdata[i].hits;
        }

        for (int i = 0; i < n_seqs; ++i) {
            total_bases += seqs[i].l_seq;
            free(seqs[i].seq);
            free(seqs[i].name);
            if (seqs[i].qual) free(seqs[i].qual);
        }
        free(seqs);
        free(threads);
        free(tdata);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double secs = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("----------------------------------------------------\n");
    printf(" Total Bases      : %lu\n", total_bases);
    printf(" Total Hits       : %lu\n", total_hits);
    printf(" Total Time       : %.4f s\n", secs);
    printf(" Throughput       : %.4f Gbp/s\n", (total_bases / 1e9) / secs);
    printf("----------------------------------------------------\n");

    mm_bseq_close(fp);
    return 0;
}
