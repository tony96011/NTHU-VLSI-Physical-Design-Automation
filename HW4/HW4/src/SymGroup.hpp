
// #pragma once
#include <string>
#include <vector>

class SymGroup {
    std::vector<std::pair<std::string,std::string>> pairs;
    std::vector<std::string> selfs;
public:
    void addSymPair(const std::string &a, const std::string &b) {
        pairs.emplace_back(a,b);
    }
    void addSymSelf(const std::string &a) {
        selfs.push_back(a);
    }
    const std::vector<std::pair<std::string,std::string>>& getSymPairs() const{
        return pairs;
    }
    const std::vector<std::string>& getSymSelfs() const{
        return selfs;
    }
    std::vector<std::string> members() const{
        auto m = selfs;
        for (auto &p : pairs) {
            m.push_back(p.first);
            m.push_back(p.second);
        }
        return m;
    }
};
