#pragma once
#include <bits/stdc++.h>
#include "BStarTree.hpp"
#include "Hardblock.hpp"
#include "SymGroup.hpp"

class Floorplan {
    int NH, NS;
    std::vector<std::string> hier;
    std::unordered_map<std::string,Hardblock> HB;
    std::unordered_map<std::string,SymGroup> SG;
    std::unordered_map<std::string,Node<int64_t>*> nodes;
    std::unordered_map<std::string,Node<int64_t>*> nodes_hier;
    std::mt19937_64 rng;            // 全局随机引擎
    std::uniform_real_distribution<double> uni01{0.0, 1.0};

    void parseInput(const std::string &);
    void writeOutput(const std::string &, int64_t, int64_t) ;
    void buildSymmetryIslands();
    void syncCoordinates(std::vector<Node<int64_t>*>& pre);
    std::pair<std::vector<Node<int64_t>*>,std::vector<Node<int64_t>*>> init_sol();
    BStarTree<int64_t> packIsland(HierNode<int64_t>* hnode);
    void perturb(std::vector<Node<int64_t>*>& pre, std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs);
    void rotate_node(std::vector<Node<int64_t>*> &pre, std::vector<std::string>& selfs);
    void checkOverlap();
    void greedyVerticalPack();
    void greedyHorizontalPack();
    void updateBoundary(int64_t &totalW, int64_t &totalH);
    void move_node(std::vector<Node<int64_t>*>& pre, std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs);
    void swap_node(std::vector<Node<int64_t>*>& pre, std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs);
    std::pair<std::vector<Node<int64_t>*>, std::vector<Node<int64_t>*>>simulatedAnnealing(const std::vector<Node<int64_t>*>& pre,const std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs);
    std::pair<std::vector<Node<int64_t>*>,std::vector<Node<int64_t>*>> clone_traversals(const std::vector<Node<int64_t>*>& pre, const std::vector<Node<int64_t>*>& in);
    void setRandomSeed(uint64_t seed);
public:
    Floorplan(const std::string &in, const std::string &out);
};