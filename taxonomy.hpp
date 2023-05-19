#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_set>

class Taxonomy
{
private:
    std::vector<int> parent;

public:
    Taxonomy() {}
    Taxonomy(char *nodes_path)
    {
        load_taxonomy(nodes_path);
    }
    void load_taxonomy(const std::string &nodes_path)
    {
        std::ifstream file(nodes_path);
        if (!file.is_open())
        {
            return;
        }
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string token;
            std::getline(iss, token, '|');
            int index = std::stoi(token);
            std::getline(iss, token, '|');
            int value = std::stoi(token);
            int min_table_size = index > value ? index : value;
            if (min_table_size >= parent.size())
            {
                parent.resize(min_table_size + 1);
            }
            parent[index] = value;
            // debug
            // if (index == 2559073) {
            //     fprintf(stderr, "son:%d, parent:%d\n", index, parent[2559073]);
            // }
            if (value == 0) {
                fprintf(stderr, "taxonomy build err: son:%d, parent:%d", index, value);
                exit(1);
            }
        }
        // fprintf(stderr, "son:%d, parent:%d\n", 2559073, parent[2559073]);
        // fprintf(stderr, "son:%d, parent:%d\n", 31957, parent[31957]);
        // exit(1);
        parent[1] = -1; // 用来做lca创建set时候的休止符
        return;
    }

    // 返回1表示所有生物，也就是没分类出来
    // 1:{10239(virus), 131567(cellular organisim), 2787823(unclassified entries), 2787854(other entries)}
    // 131567:{2(bacteria), 2157(archaea), 2759(eukaryota)}
    int lca(int p, int q)
    {
        int pp = p, qq = q;
        std::unordered_set<int> ancestors;
        while (p != -1)
        {
            if (p == q) return q; // 找到祖宗了，返回祖宗
            ancestors.insert(p);
            p = parent[p];
            if(p == 0) {
                fprintf(stderr, "lca error, parent of %d is but should not be 0, maybe a newer version of nodes.dmp can solve this problem.\n", pp);
                exit(1);
            }
        }
        while (!ancestors.count(q) && q != -1)
        {
            q = parent[q];
        }
        return q;
    }

    // lca:  如果p和q是祖孙关系，那么返回祖宗
    // lca2: 如果p和q是祖孙关系，那么返回孙子
    int lca2(int p, int q)
    {
        int pp = p, qq = q;
        std::unordered_set<int> ancestors;
        while (p != -1) // p 不断向上找祖宗
        {
            if (p == q) return pp; // 找到祖宗了，返回孙子
            ancestors.insert(p);
            p = parent[p];
            if(p == 0) {
                fprintf(stderr, "lca error, parent of %d is but should not be 0, maybe a newer version of nodes.dmp can solve this problem.\n", pp);
                exit(1);
            }
        }
        while (!ancestors.count(q) && q != -1)
        {
            q = parent[q];
        }
        return q;
    }
};
