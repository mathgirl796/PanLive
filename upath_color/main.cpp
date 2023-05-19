#include <stdio.h>
#include "utils.h"
#include "bwt.hpp"
#include "taxonomy.hpp"

char* bwt_index_path = "/home/user/duanran/HMP/ref.fa";
char* nodes_path = "/home/user/duanran/.taxonkit/nodes.dmp";
char* upath_path = "/home/user/duanran/HMP/ref.fa.ubwt.upath.fa.1";
char* output_path = "/home/user/duanran/HMP/ref.fa.ubwt.upath.fa.taxid.fa";

#define BUF_SIZ 1000000
int main() {
    fprintf(stderr, "Hello world\n");
    FILE* f = fopen(upath_path, "r");
    FILE* fo = fopen(output_path, "w");
    BWT_IDX_loader bwt_idx;
    bwt_idx.restore_idx(bwt_index_path, true, true, true, true);
    bwt_idx.bns->parse_taxid(6);
    Taxonomy taxonomy(nodes_path);
    char seq[BUF_SIZ];
    ubyte_t seq_bin[BUF_SIZ];
    uint64_t seq_num = 0;
    while (fgets(seq, BUF_SIZ, f) != NULL) {
        if(seq[0] == '>') {
            // err_puts(seq);
        }
        else {
            seq_num += 1;
            int seq_len;
            for (seq_len = 0; seq[seq_len] != '\n'; ++seq_len) {
                switch(seq[seq_len]){
                case 'A': seq_bin[seq_len] = 0; break;
                case 'C': seq_bin[seq_len] = 1; break;
                case 'G': seq_bin[seq_len] = 2; break;
                case 'T':  seq_bin[seq_len] = 3; break;
                default :
                    seq_bin[seq_len] = 4; break;
                }
            }
            seq[seq_len] = '\0';
            bwtint_t sa_begin, sa_end;
            // bwt match
            bwt_idx.bwt_match_exact_demo(seq_len, seq_bin, sa_begin, sa_end);
            // printf("sa_begin: %ld\t sa_end: %ld\n", sa_begin, sa_end);
            int final_taxid;
            for (bwtint_t i = sa_begin; i <= sa_end; ++i) {
                // fprintf(stderr, "%lu\n", i);
                bwtint_t pos = bwt_idx.bwt_sa(i);
                bool is_rev;
                bwtint_t depos = bwt_idx.sa_depos(pos, is_rev);
                uint32_t offset;
                int rid = bwt_idx.bns->bns_pos2rid(depos, offset);
                char* name = bwt_idx.bns->anns[rid].name;
                int taxid = bwt_idx.bns->anns[rid].taxid;
                if (i == sa_begin) {
                    final_taxid = taxid;
                }
                else {
                    final_taxid = taxonomy.lca(taxid, final_taxid);
                }
            }
            // fprintf(stderr, ">%d\n%s\n\n", final_taxid, seq);
            fprintf(fo, ">%d\n%s\n\n", final_taxid, seq);
            if (seq_num % 10000 == 0) {
                fprintf(stderr, "%ld\n", seq_num);
                // break;
            }
        }
    }
    fprintf(stderr, "total %ld seqs\n", seq_num);
    fclose(f);
    fclose(fo);
    fprintf(stderr, "Bye bye world\n");
}