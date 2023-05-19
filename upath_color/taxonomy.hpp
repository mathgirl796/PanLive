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
        std::vector<int> result;
        std::ifstream file(nodes_path);
        if (!file.is_open())
        {
            parent = result;
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
            if (index >= result.size())
            {
                result.resize(index + 1);
            }
            result[index] = value;
        }
        result[1] = -1; // 用来做lca创建set时候的休止符
        parent = result;
        return;
    }

    // 返回1表示所有生物，也就是没分类出来
    // 1:{10239(virus), 131567(cellular organisim), 2787823(unclassified entries), 2787854(other entries)}
    // 131567:{2(bacteria), 2157(archaea), 2759(eukaryota)}
    int lca(int p, int q)
    {
        std::unordered_set<int> ancestors;
        while (p != -1)
        {
            ancestors.insert(p);
            p = parent[p];
        }
        while (!ancestors.count(q))
        {
            q = parent[q];
        }
        return q;
    }
};
