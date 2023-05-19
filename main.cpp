#include <stdio.h>
#include "main.hpp"
#include "cmdline.hpp"

void *load_data(Shared *shared, int index)
{
    // fprintf(stderr, "load\t%d\n", index);
    int tid = index % step_num;
    Manager *m = shared->managers[tid];

    // 读取新read
    int seq_num = err_gzread(shared->fseq, m->seq, data_per_manager);
    int qual_num = err_gzread(shared->fqual, m->qual, data_per_manager);
    xassert(seq_num == qual_num, "in load data, all nums should be compat");
    // fprintf(stderr, "load index:%d\ttid:%d\t\n", index, tid);

    // 读取上个cycle的结果
    if (cycle >= 1)
    {
        xread(m->kl, sizeof(bwtint_t), seq_num * 2, shared->ftmp_last);
        xread(m->taxon, sizeof(uint32_t), seq_num, shared->ftmp_last);
        xread(m->match_len, sizeof(uint8_t), seq_num, shared->ftmp_last);
        xread(m->max_match_len, sizeof(uint8_t), seq_num, shared->ftmp_last);
    }
    else
    {
        // 第一个cycle：初始化中间结果
        for (int i = 0; i < seq_num; i++)
        {
            m->set_kl(i, 0, shared->bwt_idx->seq_len);
            m->taxon[i] = 0;
            m->match_len[i] = 0;
            m->max_match_len[i] = 0;
        }
    }
    // 判断是否为最后一块数据
    m->num = seq_num;
    if (index == 28)
    {
        fprintf(stderr, "cycle %d\tindex %d\tread %d\n", cycle, index, seq_num);
    }

    if (m->num > 0)
        return (void *)1;
    else
        return (void *)0;
}
typedef struct OccSaWorkerData
{
    Manager *m;
    int manager_index;
    int manager_tid;
} OccSaWorkerData;
void occ_sa_work(void *data, long i, int worker_id)
{    
    OccSaWorkerData *worker_data = (OccSaWorkerData *)data;
    Manager *m = worker_data->m;
    Shared *shared = (Shared *)m->shared;
    BWT_IDX_loader *bwt_idx = shared->bwt_idx;
    bntseq_t *bns = bwt_idx->bns;
    Taxonomy *taxonomy = shared->taxonomy;
    int index = worker_data->manager_index;

    bool verbose = false; if (index == 0 && i == 123) verbose = true; // debug
    if (verbose) fprintf(stderr, "index:%d, i:%ld\n", index, i);

    bool sa = false, error = false;
    bwtint_t k, l, ok, ol, sa_start, sa_end;
    uint8_t c;
    m->get_kl(i, &k, &l);
    c = 3-nst_nt4_table[m->seq[i]];

    if (verbose) { fprintf(stderr, "!before occ. k:%ld l:%ld base:%c\t", k, l, m->seq[i]); }
    // 延申一个碱基
    if (c > 3) {
        error = true;
    }
    else {
        bwt_idx->bwt_2occ(k, l, c, &ok, &ol);
        ok = bwt_idx->L2[c] + ok;
        ol = bwt_idx->L2[c] + ol;
        if (ok >= ol) {
            error = true;
        }
    }

    // 如果error，再给一次机会
    if (error) {
        bwtint_t cntk[4]; bwtint_t cntl[4];
        bwt_idx->bwt_2occ4(k, l, cntk, cntl);
        int path_N = 0;
        for(int j = 0; j < 4; j++){
            if(cntk[j] < cntl[j]){
                path_N++;
                ok = bwt_idx->L2[j] + cntk[j];
                ol = bwt_idx->L2[j] + cntl[j];
            }
        }
        if (path_N == 1) {
            error = false;
        }
    }

    if (verbose) { fprintf(stderr, "!after occ. ok:%ld ol:%ld\t", ok, ol); }

    // 判断是否需要进行sa
    if (error == true) { // 失配，可能由于N或者bwt无匹配导致
        { sa = true; sa_start = k + 1; sa_end = l + 1; }
    }
    else
    {
        // 仍能匹配
        m->match_len[i] += 1;
        m->set_kl(i, ok, ol);
        // 最后一轮执行sa
        if (cycle == read_len - 1) { sa = true; sa_start = ok + 1; sa_end = ol + 1; } 
    }

    if (verbose) { fprintf(stderr, "!bwt extend done. ok:%ld ol:%ld\n", ok, ol); }
    if (sa) {
        // 一定情况下会跳过SA
        if (sa_end - sa_start <= sa_threshold && m->match_len[i] >= match_threshold) {
            if (verbose) 
            { fprintf(stderr, "!bwt pos:%ld - %ld, match:%d, index:%d, i:%ld\n", sa_start, sa_end, m->match_len[i], index, i); }
            for (bwtint_t idx = sa_start; idx < sa_end; idx ++) {
                if (verbose) 
                { fprintf(stderr, "!!!bwt pos:%ld\t", idx); }
                bwtint_t pos = bwt_idx->bwt_sa(idx);
                bool is_rev;
                int64_t depos = bwt_idx->sa_depos(pos, is_rev);
                if (verbose) { fprintf(stderr, "!!!sadepos:done\t"); }
                uint32_t offset;
                int rid = bns->bns_pos2rid(depos, offset);
                if (verbose) { fprintf(stderr, "!!!pos2rid:done\t"); }
                int taxon = bns->anns[rid].taxid;
                if (verbose) { fprintf(stderr, "!!!old_tax:%d,new_tax:%d\t",m->taxon[i],taxon); }
                m->taxon[i] = m->taxon[i] == 0 ? taxon : taxonomy->lca2(taxon, m->taxon[i]);
                if (verbose) { fprintf(stderr, "!!!lca:done\t"); }
                if (verbose) 
                { fprintf(stderr, "!depos:%ld\tisrev:%d\ttaxon:%d\trefname:%s\n", depos, is_rev, m->taxon[i], bns->anns[rid].name); }
            }
        }
        // 无论是否执行SA，都要重置一些信息，重新开始搜索
        m->max_match_len[i] = m->max_match_len[i] > m->match_len[i] ? m->max_match_len[i] : m->match_len[i];
        m->match_len[i] = 0;
        m->set_kl(i, 0, bwt_idx->seq_len);
    }
    if (verbose)
    {
        fprintf(stderr, "match_len:%d\ttaxon:%d\t\n", m->match_len[i], m->taxon[i]);
        // exit(0);
    }
}
void occ_sa_block_work(void *data, long block_id, int worker_id)
{
    OccSaWorkerData *worker_data = (OccSaWorkerData *)data;
    long block_start = block_id * data_per_block;
    long temp_end = (block_id + 1) * data_per_block;
    long num = worker_data->m->num;
    long block_end = num < temp_end ? num : temp_end;
    for (long i = block_start; i < block_end; ++i)
    {
        occ_sa_work(data, i, worker_id);
        // fprintf(stderr, "index:%d\tblock:%ld\ti:%d\t\n", worker_data->manager_index, block_id, i);
    }
    // fprintf(stderr, "index:%d\tkt_for_worker_id:%d\tblock_s:%ld\tblock_end:%d\t\n", worker_data->manager_index, worker_id, block_start, block_end);
}
void *occ_sa(Shared *shared, int index)
{
    // fprintf(stderr, "occ\t%d\n", index);
    int tid = index % step_num;
    Manager *m = shared->managers[tid];

    // 更新kl和taxon
    OccSaWorkerData data = {m, index, tid};
    long block_num = (m->num - 1) / data_per_block + 1;
    kt_for(thread_num, occ_sa_block_work, &data, block_num);

    return (void *)1;
}
void *dump_temp(Shared *shared, int index)
{
    // fprintf(stderr, "dump\t%d\n", index);
    int tid = index % step_num;
    Manager *m = shared->managers[tid];

    // 每个cycle写出中间结果
    err_fwrite(m->kl, sizeof(bwtint_t), m->num * 2, shared->ftmp);
    err_fwrite(m->taxon, sizeof(uint32_t), m->num, shared->ftmp);
    err_fwrite(m->match_len, sizeof(uint8_t), m->num, shared->ftmp);
    err_fwrite(m->max_match_len, sizeof(uint8_t), m->num, shared->ftmp);
    // fprintf(stderr, "dump index:%d\ttid:%d\t\n", index, tid);

    if (cycle == read_len - 1)
    {
        // 最后一cycle，写出结果
        for (int i = 0; i < m->num; ++i)
        {
            fprintf(shared->foutput, "%u\n", m->taxon[i]);
        }
    }

    return (void *)1;
}
void *pipeline(void *shared, int step, void *data, int index)
{
    if (step == 0)
        return load_data((Shared *)shared, index);
    else if (step == 1)
        return occ_sa((Shared *)shared, index);
    else if (step == 2)
        return dump_temp((Shared *)shared, index);
    else
    {
        fprintf(stderr, "pipline error step:\t%d\n", step);
        return (void *)0;
    }
}

#include <string>
#include <string.h>
using namespace std;
int main(int argc, char **argv)
{
    double start_cpu_time = cputime(), start_real_time = realtime();
    double __cpu_time, __real_time;
    __cpu_time = cputime(); __real_time = realtime();

    cmdline::parser a;
    a.add<int>("thread_num", 'j', "thread num", false, 4);
    a.add<int>("data_per_manager", 'm', "data per manager", false, 5000);
    a.add<int>("data_per_block", 'b', "data per block in a manager", false, 1000);
    a.add<int>("read_len", 'l', "length of read", true);
    a.add<int>("sa_threshold", 0, "max sa range when need to do sa", false, 1000000);
    a.add<int>("match_threshold", 0, "min match num when need to do sa", false, 0);
    a.add<string>("input_dir", 'i', "input SAMPLE directory", true);
    a.add<string>("idx_path", 'x', "path to bwt index prefix", true);
    a.add<string>("tmp_dir", 't', "path to tmp dir to store temp files", false, "./");
    a.add<string>("output_path", 'o', "output file path", true);
    a.add<string>("nodes_path", 'n', "path to nodes.dmp", true);
    a.parse_check(argc, argv);

    thread_num = a.get<int>("thread_num");
    data_per_manager = a.get<int>("data_per_manager");
    data_per_block = a.get<int>("data_per_block");
    read_len = a.get<int>("read_len");
    sa_threshold = a.get<int>("sa_threshold");
    match_threshold = a.get<int>("match_threshold");
    input_dir = strdup(a.get<string>("input_dir").c_str());
    idx_path = strdup(a.get<string>("idx_path").c_str());
    tmp_dir = strdup(a.get<string>("tmp_dir").c_str());
    output_path = strdup(a.get<string>("output_path").c_str());
    nodes_path = strdup(a.get<string>("nodes_path").c_str());

    fprintf(stderr, "Hello world!\n");
    Shared shared;
    shared.init();
    char fn[1024];

    // open result file to write final result
    sprintf(fn, "%s", output_path);
    shared.foutput = xopen(fn, "wb");

    fprintf(stderr, "init over: CPU time: %.3f sec; real time: %.3f sec\n", cputime() - __cpu_time, realtime() - __real_time);
    for (cycle = 0; cycle < read_len; cycle++)
    {
        __cpu_time = cputime(); __real_time = realtime();
        char log_title[1024]; 
        sprintf(log_title, "[cycle %d]", cycle);
        fprintf(stderr, "%s start.\n", log_title);
        // open new read and new qual
        sprintf(fn, "%s/S_%03d.bin.gz", input_dir, cycle);
        shared.fseq = xzopen(fn, "rb");
        sprintf(fn, "%s/Q_%03d.bin.gz", input_dir, cycle);
        shared.fqual = xzopen(fn, "rb");

        // open new tmp to write tmp result
        sprintf(fn, "%s/%03d.tmp", tmp_dir, cycle);
        shared.ftmp = xopen(fn, "wb");
        if (cycle >= 1)
        {
            // open last temp to read tmp result
            sprintf(fn, "%s/%03d.tmp", tmp_dir, cycle - 1);
            shared.ftmp_last = xopen(fn, "rb");
        }

        // run pipeline for this cycle
        kt_pipeline(step_num, pipeline, &shared, step_num);

        // release resources
        err_gzclose(shared.fseq);
        err_gzclose(shared.fqual);
        err_fclose(shared.ftmp);
        if (cycle >= 1)
        {
            err_fclose(shared.ftmp_last);
            // delete last cycle tmp file
            sprintf(fn, "%s/%03d.tmp", tmp_dir, cycle - 1);
            // fprintf(stderr, "[main] remove %s\n", fn);
            xrm(fn);
        }
        fprintf(stderr, "%s end. CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
    }
    // delete final cycle tmp file
    sprintf(fn, "%s/%03d.tmp", tmp_dir, cycle - 1);
    xrm(fn);

    // release resource
    err_fclose(shared.foutput);

    fprintf(stderr, "[ALL] end. total CPU time: %.3f sec; real time: %.3f sec\n", cputime() - start_cpu_time, realtime() - start_real_time);

    // stop program
    fprintf(stderr, "Bye bye world!\n");
    return 0;
}