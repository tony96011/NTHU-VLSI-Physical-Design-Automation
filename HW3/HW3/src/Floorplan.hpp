#include <bits/stdc++.h>

#ifndef NET_HPP
#define NET_HPP
#include "Net.hpp"
#endif

#ifndef HARDBLOCK_HPP
#define HARDBLOCK_HPP
#include "Hardblock.hpp"
#endif

#ifndef PIN_HPP
#define PIN_HPP
#include "Pin.hpp"
#endif

#ifndef COORD_HPP
#define COORD_HPP
#include "Coord.hpp"
#endif

#ifndef NODE_HPP
#define NODE_HPP
#include "Node.hpp"
#endif


class Floorplan {
public:
    // Variables
    unsigned int seed;
    int num_hardblocks;
    int num_terminals;
    int num_nets;
    int total_area;
    long long best_w, best_h;
    double dead_space_ratio;
    long long max_coord_x, max_coord_y;
    std::unordered_map<std::string, Hardblock> hardblocks;
    std::unordered_map<std::string, Pin> pins;
    std::vector<Net> nets;
    std::chrono::high_resolution_clock::time_point start;

    // Constructors
    Floorplan();
    Floorplan(std::string input_file, std::string output, double ratio);

    // Functions
    void set_seed(std::string input_file);
    void read_input(std::string input_file);
    void write_floorplan(std::string filename);
    void calc_max_coord();
    void calc_total_area();
    void update_coord(std::vector<std::vector<Node>>& record, int index, int min_at);
    void swap_adjacent_operand(std::vector<std::string>& sol);
    void invert_chain(std::vector<std::string>& sol);
    void swap_operand_operator(std::vector<std::string>& sol);
    void swap_random_operand(std::vector<std::string>& sol);
    int get_wirelength();
    long long get_cost(std::vector<std::string> sol, bool only_area);
    // int get_cost(std::vector<std::string> sol, bool only_area);
    int get_area(std::vector<std::string>& sol, long long& w, long long& h);
    std::vector<std::string> init_sol();
    std::vector<std::string> gen_neighbor(std::vector<std::string> sol, int r);
    double calc_initial_temp(std::vector<std::string>& sol, double p, bool only_area);
    std::vector<std::string> simulated_annealing(std::vector<std::string>& sol, bool only_area);
    std::vector<Node> stockmeyer(std::vector<Node> l, std::vector<Node> r, std::string type, int index);
};