#include <bits/stdc++.h>
#include "Floorplan.hpp"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] 
                  << " <input_file.txt> <output_file.out> <dead_space_ratio>\n";
        return 1;
    }
    //計算時間
    auto start = std::chrono::high_resolution_clock::now();
    Floorplan floorplan(argv[1], argv[2], std::stod(argv[3]));
    
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}
