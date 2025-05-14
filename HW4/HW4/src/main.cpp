#include <bits/stdc++.h>
#include "Floorplan.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr<<"Usage: "<<argv[0]<<" <input> <output>\n";
        return 1;
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    Floorplan fp(argv[1], argv[2]);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    std::cout<<"Elapsed: "<<dt.count()<<" s\n";
    return 0;
}