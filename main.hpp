#pragma once

#include <stdint.h>
#include "utils.h"
#include "cmdline.hpp"
#include "bwt.hpp"
#include "kthread.h"
#include "taxonomy.hpp"

// user parameters
int thread_num;
int data_per_manager;
int data_per_block;
int read_len;
int sa_threshold;
int match_threshold;
char *input_dir;
char *idx_path;
char *tmp_dir;
char *output_path;
char *nodes_path;

// program used global vals
int step_num = 3;
int cycle;

typedef struct Manager
{
    void *shared;
    uint8_t *seq;
    uint8_t *qual;
    bwtint_t *kl;
    uint32_t *taxon;
    uint8_t *match_len;
    uint8_t *max_match_len;
    int num;
    void init()
    {
        this->seq = (uint8_t *)xcalloc(data_per_manager, sizeof(uint8_t));
        this->qual = (uint8_t *)xcalloc(data_per_manager, sizeof(uint8_t));
        this->kl = (bwtint_t *)xcalloc(data_per_manager * 2, sizeof(bwtint_t));
        this->taxon = (uint32_t *)xcalloc(data_per_manager, sizeof(int));
        this->match_len = (uint8_t *)xcalloc(data_per_manager, sizeof(uint8_t));
        this->max_match_len = (uint8_t *)xcalloc(data_per_manager, sizeof(uint8_t));
    }
    void set_kl(int i, bwtint_t k, bwtint_t l)
    {
        this->kl[i * 2] = k;
        this->kl[i * 2 + 1] = l;
    }
    void get_kl(int i, bwtint_t *k, bwtint_t *l)
    {
        *k = this->kl[i * 2];
        *l = this->kl[i * 2 + 1];
    }
} Manager;

typedef struct Shared
{
    Manager **managers;
    BWT_IDX_loader *bwt_idx;
    Taxonomy *taxonomy;
    gzFile fseq;
    gzFile fqual;
    FILE *ftmp;
    FILE *ftmp_last;
    FILE *foutput;
    void init()
    {
        this->managers = (Manager **)xcalloc(step_num, sizeof(Manager *));
        this->bwt_idx = (BWT_IDX_loader *)xcalloc(1, sizeof(BWT_IDX_loader));
        this->bwt_idx->restore_idx(idx_path, true, true, true, true);
        this->bwt_idx->bns->parse_taxid(0);
        this->taxonomy = (Taxonomy *)xcalloc(1, sizeof(Taxonomy));
        this->taxonomy->load_taxonomy(nodes_path);
        for (int i = 0; i < step_num; ++i)
        {
            this->managers[i] = (Manager *)xcalloc(1, sizeof(Manager));
            this->managers[i]->init();
            this->managers[i]->shared = this;
        }
    }
} Shared;

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};