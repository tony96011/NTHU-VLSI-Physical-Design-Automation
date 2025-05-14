#include <bits/stdc++.h>
#include <chrono>
using namespace std;
using namespace chrono;

class LibCell {
public:
    string name;
    long width, height;
    LibCell(string n = "", long w = 0, long h = 0) : name(n), width(w), height(h) {}
};

class Tech {
public:
    string name;
    unordered_map<string, LibCell> libCells; // key: LibCell name, value: LibCell object
    Tech(string n = "") : name(n) {}
};

class Cell {
public:
    string name;
    string libCell;
    int techIndex;         // Index into the vector<Tech>
    bool isA = false, isLocked = false;
    vector<int> netIndices;  // Indices into the vector<Net>
    int gain = 0;
    // cluster-related members omitted for simplicity
    // originalCells can be used to record indices if needed

    Cell(string n = "", string lc = "", int tIndex = -1, bool a = false, bool lock = false)
        : name(n), libCell(lc), techIndex(tIndex), isA(a), isLocked(lock) {}

    double area(const vector<Tech>& techs) const {
        const Tech &tech = techs[techIndex];
        auto it = tech.libCells.find(libCell);
        return (it != tech.libCells.end()) 
            ? static_cast<double>(it->second.width * it->second.height)
            : 0.0;
    }
};

class Net {
public:
    string name;
    vector<int> cellIndices; // Indices into the vector<Cell>
    int weight;
    long numA = 0, numB = 0;
    Net(string n = "", int w = 0) : name(n), weight(w) {}
};

class Die {
public:
    string name;
    int techIndex; // Index into vector<Tech>
    long width, height;
    unordered_map<string, int> cellMap; // key: cell name, value: index in vector<Cell>
    float utilization; // in percent

    Die(string n = "", int tIndex = -1, long w = 0, long h = 0, float u = 0.0f)
        : name(n), techIndex(tIndex), width(w), height(h), utilization(u) {}

    double maxArea() const { return width * height * (utilization / 100.0); }
    double totalArea() const { return width * height; }
    double currentArea(const vector<Cell>& cells, const vector<Tech>& techs) const {
        double sum = 0.0;
        for (const auto& kv : cellMap)
            sum += cells.at(kv.second).area(techs);
        return sum;
    }
    bool withinUtilization(const vector<Cell>& cells, const vector<Tech>& techs) const { 
        return currentArea(cells, techs) <= maxArea(); 
    }
};

//
// processInput: Reads the input file and builds the Tech, LibCell, Cell, and Net objects,
// setting up the relationships among them. Notice that all relationships are stored as indices
// into the corresponding vectors rather than using raw pointers.
//
float getCutSize(const vector<Net>& nets, const vector<Cell>& cells) {
    float cutSize = 0.0f;
    for (const auto& net : nets) {
        bool hasA = false, hasB = false;
        for (int cellIndex : net.cellIndices) {
            if (cells[cellIndex].isA)
                hasA = true;
            else
                hasB = true;
            if (hasA && hasB) {
                cutSize += net.weight;
                break;
            }
        }
    }
    return cutSize;
}

void printDieCells(const Die& dieA, const Die& dieB) {
    cout << "Cells in " << dieA.name << ":\n";
    for (const auto& kv : dieA.cellMap) {
        cout << kv.first << " ";
    }
    cout << "\n\nCells in " << dieB.name << ":\n";
    for (const auto& kv : dieB.cellMap) {
        cout << kv.first << " ";
    }
    cout << "\n";
}

void printGainBuckets(const map<int, unordered_set<int>, greater<int>>& gainA,
    const map<int, unordered_set<int>, greater<int>>& gainB,
    const vector<Cell>& cells) {
    cout << "GainA Buckets:\n";
    for (const auto& bucket : gainA) {
    cout << "Gain " << bucket.first << ": ";
        for (int cellIdx : bucket.second) {
            cout << cells[cellIdx].name << " ";
        }
        cout << "\n";
    }

    cout << "\nGainB Buckets:\n";
    for (const auto& bucket : gainB) {
        cout << "Gain " << bucket.first << ": ";
        for (int cellIdx : bucket.second) {
            cout << cells[cellIdx].name << " ";
        }
        cout << "\n";
    }
}

pair<double, double> processInput(const string& file,
                                  unordered_map<string, int>& techMap, 
                                  vector<Tech>& techs,
                                  Die& dieA, Die& dieB,
                                  vector<Cell>& cells, 
                                  vector<Net>& nets) {
    ifstream in(file);
    if (!in) throw runtime_error("Cannot open input file: " + file);

    double totalAreaA = 0.0, totalAreaB = 0.0;
    long width = 0, height = 0, numCells = 0, numNets = 0;
    unordered_map<string, int> cellMap;  // maps cell name to index in vector<Cell>

    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string keyword;
        ss >> keyword;
        int currentTechIndex = -1;  // Index of the current Tech object being processed

        if (keyword == "Tech") {
            int num_libcells;
            string name; 
            ss >> name >> num_libcells;
            techs.push_back(Tech(name));
            currentTechIndex = techs.size() - 1;  // update current tech index
            techMap[name] = currentTechIndex;
            // Read each LibCell from the next lines
            for (int i = 0; i < num_libcells; i++) {
                if(!getline(in, line)) break; // read next LibCell line
                stringstream libSS(line);
                string libcell_keyword, libcell_name;
                int w, h;
                libSS >> libcell_keyword >> libcell_name >> w >> h;
                techs[currentTechIndex].libCells[libcell_name] = LibCell(libcell_name, w, h);
            }
        }
        else if (keyword == "DieSize") {
            ss >> width >> height;
        }
        else if (keyword == "DieA") {
            string techName; int util;
            ss >> techName >> util;
            int tIndex = techMap[techName];
            dieA = Die("DieA", tIndex, width, height, util);
        }
        else if (keyword == "DieB") {
            string techName; int util;
            ss >> techName >> util;
            int tIndex = techMap[techName];
            dieB = Die("DieB", tIndex, width, height, util);
        }
        else if (keyword == "NumCells") {
            ss >> numCells;
        }
        else if (keyword == "Cell" && numCells > 0) {
            string name, libCell;
            ss >> name >> libCell;
            int initTech = dieA.techIndex;
            Cell cell(name, libCell, initTech, true, false);
            // Calculate area in DieA and DieB (by temporarily switching tech index)
            cell.techIndex = dieA.techIndex;
            double areaA = cell.area(techs);
            cell.techIndex = dieB.techIndex;
            double areaB = cell.area(techs);
            cell.isA = (areaA <= areaB);
            cell.techIndex = cell.isA ? dieA.techIndex : dieB.techIndex;
            cells.push_back(cell);
            int idx = cells.size() - 1;
            cellMap[name] = idx;
            if (cell.isA) {
                dieA.cellMap[name] = idx;
                totalAreaA += cells[idx].area(techs);
            } else {
                dieB.cellMap[name] = idx;
                totalAreaB += cells[idx].area(techs);
            }
            numCells--;
        }
        else if (keyword == "NumNets") {
            ss >> numNets;
        }
        else if (keyword == "Net" && numNets > 0) {
            string name; int weight; long cellCount;
            ss >> name >> cellCount >> weight;
            nets.push_back(Net(name, weight));
            int netIdx = nets.size() - 1;
            for (long i = 0; i < cellCount && getline(in, line); i++) {
                stringstream cellSS(line);
                string cellKeyword, cellName;
                cellSS >> cellKeyword >> cellName;
                if (cellKeyword != "Cell" || cellMap.find(cellName) == cellMap.end())
                    continue;
                int cIdx = cellMap[cellName];
                nets[netIdx].cellIndices.push_back(cIdx);
                // Add the net index to the corresponding cell
                cells[cIdx].netIndices.push_back(netIdx);
                if (cells[cIdx].isA)
                    nets[netIdx].numA++;
                else
                    nets[netIdx].numB++;
            }
            numNets--;
        }
    }

    cout << "Input processed: " << cells.size() << " cells, " << nets.size() << " nets\n";
    cout << "DieA utilization: " << (totalAreaA / dieA.totalArea() * 100.0) << "%\n";
    cout << "DieB utilization: " << (totalAreaB / dieB.totalArea() * 100.0) << "%\n";
    return {totalAreaA, totalAreaB};
}

// Cluster function using index-based data structures.
void cluster(vector<Cell>& cells, const vector<Tech>& techs,
             Die& dieA, Die& dieB, vector<Net>& nets,
             double& totalAreaA, double& totalAreaB) {
    cout << "Running Clustering...\n";
    unordered_set<int> visited;  // indices of cells that are clustered
    dieA.cellMap.clear();
    dieB.cellMap.clear();

    // Create a list of net indices sorted by descending net weight.
    vector<int> sortedNetIndices;
    for (int i = 0; i < nets.size(); i++) {
        sortedNetIndices.push_back(i);
    }
    sort(sortedNetIndices.begin(), sortedNetIndices.end(), [&](int i, int j) {
        return nets[i].weight > nets[j].weight;
    });

    // Define a Cluster structure that holds indices of cells plus the computed areas
    struct Cluster {
        vector<int> cellIndices;
        double areaA, areaB; // total area if cluster is assigned to DieA and DieB respectively
        bool preferA() const { return areaA < areaB; }
        double areaDiff() const { return fabs(areaA - areaB); }
    };
    vector<Cluster> clusters;

    // Helper lambda: computes how much the area would change if a cell were moved.
    auto areaDiffLambda = [&](int cellIndex) -> double {
        double currentArea = cells[cellIndex].area(techs);
        int currentTech = cells[cellIndex].techIndex;
        // simulate moving: if currently in A then compute area in B; else in A.
        int newTech = cells[cellIndex].isA ? dieB.techIndex : dieA.techIndex;
        cells[cellIndex].techIndex = newTech;
        double newArea = cells[cellIndex].area(techs);
        cells[cellIndex].techIndex = currentTech;
        return newArea - currentArea;
    };

    // Form clusters by iterating over nets (starting with the highest weight).
    for (int netIdx : sortedNetIndices) {
        const Net &net = nets[netIdx];
        vector<int> candidateIndices;
        for (int cellIdx : net.cellIndices) {
            if (visited.find(cellIdx) == visited.end())
                candidateIndices.push_back(cellIdx);
        }
        // Only form a cluster if all connected cells are unvisited.
        if (candidateIndices.size() < net.cellIndices.size())
            continue;
        double areaA = 0.0, areaB = 0.0;
        for (int cellIdx : candidateIndices) {
            visited.insert(cellIdx);
            // Calculate area assuming assignment to DieA.
            cells[cellIdx].techIndex = dieA.techIndex;
            areaA += cells[cellIdx].area(techs);
            // Calculate area assuming assignment to DieB.
            cells[cellIdx].techIndex = dieB.techIndex;
            areaB += cells[cellIdx].area(techs);
            // Restore current tech assignment.
            cells[cellIdx].techIndex = cells[cellIdx].isA ? dieA.techIndex : dieB.techIndex;
        }
        clusters.push_back({candidateIndices, areaA, areaB});
    }

    // For any cell that hasn't been clustered, create a single-cell cluster.
    for (int i = 0; i < cells.size(); i++) {
        if (visited.find(i) != visited.end())
            continue;
        cells[i].techIndex = dieA.techIndex;
        double areaA = cells[i].area(techs);
        cells[i].techIndex = dieB.techIndex;
        double areaB = cells[i].area(techs);
        cells[i].techIndex = cells[i].isA ? dieA.techIndex : dieB.techIndex;
        clusters.push_back({vector<int>{i}, areaA, areaB});
        visited.insert(i);
    }

    // Sort clusters in descending order by the absolute area difference.
    sort(clusters.begin(), clusters.end(), [&](const Cluster &a, const Cluster &b) {
        return a.areaDiff() > b.areaDiff();
    });

    // Reset total areas
    totalAreaA = totalAreaB = 0.0;
    // For each cluster, assign all its cells to the preferred die if possible.
    for (auto &cluster : clusters) {
        bool toA = cluster.preferA(); // assign to DieA if areaA < areaB
        double area = toA ? cluster.areaA : cluster.areaB;
        Die &targetDie = toA ? dieA : dieB;
        double &targetArea = toA ? totalAreaA : totalAreaB;
        double &altArea = toA ? totalAreaB : totalAreaA;
        double alt = toA ? cluster.areaB : cluster.areaA;

        // Lambda to assign a cluster's cells to a given die.
        auto assignCluster = [&](Die &die, bool assignToA, double &areaSum) {
            for (int cellIdx : cluster.cellIndices) {
                cells[cellIdx].isA = assignToA;
                cells[cellIdx].techIndex = die.techIndex;
                die.cellMap[cells[cellIdx].name] = cellIdx;
                areaSum += cells[cellIdx].area(techs);
            }
        };

        if (targetArea + area <= targetDie.maxArea())
            assignCluster(targetDie, toA, targetArea);
        else if (altArea + alt <= (toA ? dieB : dieA).maxArea())
            assignCluster(toA ? dieB : dieA, !toA, altArea);
        else {
            // As a fallback, choose the die with the smaller computed area.
            if (cluster.areaA < cluster.areaB)
                assignCluster(dieA, true, totalAreaA);
            else
                assignCluster(dieB, false, totalAreaB);
        }
        // Optionally, lock one cell in the cluster to preserve its assignment.
        if (!cluster.cellIndices.empty())
            cells[cluster.cellIndices[0]].isLocked = true;
    }

    // Adjust utilization if either die is over the allowed limit.
    auto adjustUtilization = [&](Die &from, Die &to, double &fromArea, double &toArea) {
        vector<int> cellIndices;
        for (const auto &kv : from.cellMap)
            cellIndices.push_back(kv.second);
        sort(cellIndices.begin(), cellIndices.end(), [&](int a, int b) {
            return areaDiffLambda(a) < areaDiffLambda(b);
        });
        double maxUtilRatio = from.utilization / 100.0;
        for (int cellIdx : cellIndices) {
            if (fromArea / from.totalArea() <= maxUtilRatio)
                break;
            from.cellMap.erase(cells[cellIdx].name);
            cells[cellIdx].isA = !cells[cellIdx].isA;
            fromArea -= cells[cellIdx].area(techs);
            cells[cellIdx].techIndex = to.techIndex;
            to.cellMap[cells[cellIdx].name] = cellIdx;
            toArea += cells[cellIdx].area(techs);
            // Update net counts as needed (omitted here for brevity).
        }
    };

    // Iteratively adjust until both dies meet the utilization constraints.
    while (totalAreaA / dieA.totalArea() > dieA.utilization / 100.0 ||
           totalAreaB / dieB.totalArea() > dieB.utilization / 100.0) {
        if (totalAreaA / dieA.totalArea() > dieA.utilization / 100.0)
            adjustUtilization(dieA, dieB, totalAreaA, totalAreaB);
        if (totalAreaB / dieB.totalArea() > dieB.utilization / 100.0)
            adjustUtilization(dieB, dieA, totalAreaB, totalAreaA);
    }

    // Finally, update each net’s counts based on the final cell assignments.
    for (auto &net : nets) {
        net.numA = net.numB = 0;
        for (int cellIdx : net.cellIndices) {
            if (cells[cellIdx].isA)
                net.numA++;
            else
                net.numB++;
        }
    }
}

int computeGain(vector<Cell>& cells, const vector<Net>& nets, map<int, unordered_set<int>, greater<int>>& gainA, map<int, unordered_set<int>, greater<int>>& gainB) {
    cout << "Computing gains...\n";
    int maxGain = INT_MIN;
    gainA.clear();
    gainB.clear();

    for (int i = 0; i < cells.size(); i++) {
        if (cells[i].isLocked) continue;
        int gain = 0;
        // Iterate through all nets connected to the cell.
        for (int netIndex : cells[i].netIndices) {
            const Net &net = nets[netIndex];
            if (cells[i].isA) {
                if (net.numA == 1) gain += net.weight;
                if (net.numB == 0) gain -= net.weight;
            } else {
                if (net.numB == 1) gain += net.weight;
                if (net.numA == 0) gain -= net.weight;
            }
        }
        cells[i].gain = gain;
        maxGain = max(maxGain, gain);
        // Insert the cell's index into the appropriate gain bucket.
        if (cells[i].isA)
            gainA[gain].insert(i);
        else
            gainB[gain].insert(i);
    }
    return maxGain;
}


// Returns the index of the cell with maximum gain from the target gain bucket.
// If a cell is selected, it is locked.
int getMaxGainCell(const vector<Cell>& cells, const vector<Tech>& techs,
                   const map<int, unordered_set<int>, greater<int>>& gainA,
                   const map<int, unordered_set<int>, greater<int>>& gainB,
                   const Die& dieA, const Die& dieB,
                   double totalAreaA, double totalAreaB, int maxGain) {
    int maxGainA = gainA.empty() ? INT_MIN : gainA.begin()->first;
    int maxGainB = gainB.empty() ? INT_MIN : gainB.begin()->first;
    bool tie = (maxGainA == maxGainB && maxGainA == maxGain);
    bool useA = tie ? ((dieA.maxArea() - totalAreaA) < (dieB.maxArea() - totalAreaB))
                    : (maxGainA > maxGainB);

    const auto &targetGain = useA ? gainA : gainB;
    int bestCellIdx = -1;
    double minArea = numeric_limits<double>::max();

    auto it = targetGain.find(maxGain);
    if (it != targetGain.end()) {
        for (int idx : it->second) {
            if (cells[idx].isLocked) continue;
            int origTech = cells[idx].techIndex;
            int newTech = cells[idx].isA ? dieB.techIndex : dieA.techIndex;
            // Create a temporary copy to simulate the move.
            Cell temp = cells[idx];
            temp.techIndex = newTech;
            double area = temp.area(techs);
            if (area < minArea) {
                minArea = area;
                bestCellIdx = idx;
            }
        }
    }
    // Lock the best candidate.
    if (bestCellIdx != -1) {
        // Note: Since cells is passed by const reference here, in fm we
        // will update the actual cell once we have its index.
        // For now, we return the index so that the caller can lock it.
    }
    return bestCellIdx;
}

// Attempts to move the cell with index cellIdx from one die to the other.
// Updates the corresponding gain bucket maps and total area variables.
// Returns true if the cell was successfully moved.
bool moveCell(int cellIdx, vector<Cell>& cells, const vector<Tech>& techs,
              Die& dieA, Die& dieB,
              map<int, unordered_set<int>, greater<int>>& gainA,
              map<int, unordered_set<int>, greater<int>>& gainB,
              bool lock, double& totalAreaA, double& totalAreaB, int& maxGain) {
    if (cellIdx < 0) return false;
    Cell &cell = cells[cellIdx];
    bool fromA = cell.isA;
    Die &fromDie = fromA ? dieA : dieB;
    Die &toDie   = fromA ? dieB : dieA;
    double &fromArea = fromA ? totalAreaA : totalAreaB;
    double &toArea   = fromA ? totalAreaB : totalAreaA;

    double oldArea = cell.area(techs);
    // Simulate moving by switching tech assignment.
    cell.techIndex = toDie.techIndex;
    double newArea = cell.area(techs);

    if (toArea + newArea > toDie.maxArea()) {
        cell.techIndex = fromDie.techIndex; // revert
        cell.isLocked = true;
        auto &fromGain = fromA ? gainA : gainB;
        fromGain[cell.gain].erase(cellIdx);
        if (fromGain[cell.gain].empty()) fromGain.erase(cell.gain);
        maxGain = max(gainA.empty() ? INT_MIN : gainA.begin()->first,
                      gainB.empty() ? INT_MIN : gainB.begin()->first);
        return false;
    }

    // Remove cell from its current gain bucket.
    auto &fromGain = fromA ? gainA : gainB;
    fromGain[cell.gain].erase(cellIdx);
    if (fromGain[cell.gain].empty()) fromGain.erase(cell.gain);

    // Update die assignments.
    fromDie.cellMap.erase(cell.name);
    toDie.cellMap[cell.name] = cellIdx;
    fromArea -= oldArea;
    toArea += newArea;
    // Flip partition.
    cell.isA = !cell.isA;
    cell.isLocked = lock;

    maxGain = max(gainA.empty() ? INT_MIN : gainA.begin()->first,
                  gainB.empty() ? INT_MIN : gainB.begin()->first);
    return true;
}

void updateGain(int cellIdx, vector<Cell>& cells, vector<Net>& nets,
                map<int, unordered_set<int>, greater<int>>& gainA,
                map<int, unordered_set<int>, greater<int>>& gainB,
                int& maxGain, const vector<Tech>& techs) {
    if (cellIdx == -1) return;
    Cell& cell = cells[cellIdx];
    bool fromA = !cell.isA; // Original partition before move
    int oldGain = cell.gain;

    // Update net counters for all nets connected to this cell (sequentially)
    for (int netIdx : cell.netIndices) {
        Net& net = nets[netIdx];
        if (fromA) {
            net.numA--;
            net.numB++;
        } else {
            net.numA++;
            net.numB--;
        }
    }
    cell.gain = -oldGain;

    // Parallel region for updating gains of neighboring cells.
    #pragma omp parallel
    {
        // Each thread keeps its own local gain maps and local maxGain.
        map<int, unordered_set<int>, greater<int>> localGainA;
        map<int, unordered_set<int>, greater<int>> localGainB;
        int localMaxGain = INT_MIN;

        // Parallelize over the nets connected to the moved cell.
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < cell.netIndices.size(); ++i) {
            int netIdx = cell.netIndices[i];
            Net& net = nets[netIdx];
            long F_new = fromA ? net.numA : net.numB;
            long T = fromA ? net.numB - 1 : net.numA - 1;

            // For each cell in the net, update its gain.
            for (int cIdx : net.cellIndices) {
                if (cIdx == cellIdx || cells[cIdx].isLocked)
                    continue;
                int oldGainC;
                // Read current gain in a critical section.
                #pragma omp critical(read_gain)
                {
                    oldGainC = cells[cIdx].gain;
                }
                int deltaGain = 0;
                if (fromA) {
                    if (cells[cIdx].isA && T == 0) deltaGain += net.weight;
                    else if (!cells[cIdx].isA && T == 1) deltaGain -= net.weight;
                    if (!cells[cIdx].isA && F_new == 0) deltaGain -= net.weight;
                    else if (cells[cIdx].isA && F_new == 1) deltaGain += net.weight;
                } else {
                    if (!cells[cIdx].isA && T == 0) deltaGain += net.weight;
                    else if (cells[cIdx].isA && T == 1) deltaGain -= net.weight;
                    if (cells[cIdx].isA && F_new == 0) deltaGain -= net.weight;
                    else if (!cells[cIdx].isA && F_new == 1) deltaGain += net.weight;
                }

                if (deltaGain != 0) {
                    // Update the gain for cell cIdx within a critical section.
                    #pragma omp critical(update_gain)
                    {
                        int currentGain = cells[cIdx].gain;
                        cells[cIdx].gain = currentGain + deltaGain;
                        // Update the local gain map for this thread.
                        auto& localMap = cells[cIdx].isA ? localGainA : localGainB;
                        if (localMap.find(currentGain) != localMap.end()) {
                            localMap[currentGain].erase(cIdx);
                            if (localMap[currentGain].empty())
                                localMap.erase(currentGain);
                        }
                        localMap[cells[cIdx].gain].insert(cIdx);
                        localMaxGain = max(localMaxGain, cells[cIdx].gain);
                    }
                }
            }
        }

        // Merge local gain maps and local maxGain into global maps.
        #pragma omp critical(merge_maps)
        {
            for (const auto& [gain, cellSet] : localGainA) {
                gainA[gain].insert(cellSet.begin(), cellSet.end());
            }
            for (const auto& [gain, cellSet] : localGainB) {
                gainB[gain].insert(cellSet.begin(), cellSet.end());
            }
            maxGain = max(maxGain, localMaxGain);
        }
    }

    // Final computation of maxGain based on global gain maps.
    maxGain = max(gainA.empty() ? INT_MIN : gainA.begin()->first,
                  gainB.empty() ? INT_MIN : gainB.begin()->first);
}


// FM algorithm: Iteratively move cells to reduce the cut size.
// Returns the maximum sum of gains achieved.
int fm(vector<Cell>& cells, vector<Net>& nets, const vector<Tech>& techs,
       Die& dieA, Die& dieB,
       map<int, unordered_set<int>, greater<int>>& gainA,
       map<int, unordered_set<int>, greater<int>>& gainB,
       double& totalAreaA, double& totalAreaB, int& maxGain) {
    cout << "Running FM...\n";
    int currentSum = 0, maxSum = 0, maxStep = 0;
    stack<int> moved;
    int cellIdx = getMaxGainCell(cells, techs, gainA, gainB, dieA, dieB, totalAreaA, totalAreaB, maxGain);
    long count = cells.size();

    while (cellIdx != -1 && count--) {
        if (moveCell(cellIdx, cells, techs, dieA, dieB, gainA, gainB, true, totalAreaA, totalAreaB, maxGain)) {
            moved.push(cellIdx);
            currentSum += cells[cellIdx].gain;
            updateGain(cellIdx, cells, nets, gainA, gainB, maxGain, techs);
            if (currentSum > maxSum) {
                maxSum = currentSum;
                maxStep = moved.size();
            }
        }
        cellIdx = getMaxGainCell(cells, techs, gainA, gainB, dieA, dieB, totalAreaA, totalAreaB, maxGain);
        if (currentSum < 0 || cells[cellIdx].gain < 0) 
            break;
    }

    // Roll back moves beyond the best partial sum.
    for (long i = moved.size() - maxStep; i > 0; i--) {
        cellIdx = moved.top();
        moved.pop();
        moveCell(cellIdx, cells, techs, dieA, dieB, gainA, gainB, false, totalAreaA, totalAreaB, maxGain);
        updateGain(cellIdx, cells, nets, gainA, gainB, maxGain, techs);
    }
    return maxSum;
}



void processOutput(const string& file, const Die& dieA, const Die& dieB, int cutSize) {
    ofstream out(file);
    if (!out)
        throw runtime_error("Cannot open output file: " + file);

    out << "CutSize " << cutSize << "\n";
    out << "DieA " << dieA.cellMap.size() << "\n";
    for (const auto& pair : dieA.cellMap)
        out << pair.first << "\n";
    out << "DieB " << dieB.cellMap.size() << "\n";
    for (const auto& pair : dieB.cellMap)
        out << pair.first << "\n";

    cout << "Output written to " << file << ".\n";
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    auto start = high_resolution_clock::now();
    string inputFile = argv[1], outputFile = argv[2];

    unordered_map<string, int> techMap;  // maps Tech name to index in vector<Tech>
    vector<Tech> techs;
    Die dieA, dieB;
    vector<Cell> cells;
    vector<Net> nets;

    cout << "Reading input...\n";
    auto areas = processInput(inputFile, techMap, techs, dieA, dieB, cells, nets);
    // cutsize
    // printDieCells(dieA, dieB);
    float cutSize = getCutSize(nets, cells);
    cout << "CutSize: " << cutSize << "\n";

    double totalAreaA = areas.first, totalAreaB = areas.second;
    cluster(cells, techs, dieA, dieB, nets, totalAreaA, totalAreaB);

    map<int, unordered_set<int>, greater<int>> gainA, gainB;
    int maxGain = computeGain(cells, nets, gainA, gainB);
    
    // print cells in dieA and dieB
    // printDieCells(dieA, dieB);
    // printGainBuckets(gainA, gainB, cells);

    int totalGain;
    int pass = 2;
    do {
        totalGain = fm(cells, nets, techs, dieA, dieB, gainA, gainB, totalAreaA, totalAreaB, maxGain);
        for (auto& cell : cells) cell.isLocked = false;
        maxGain = computeGain(cells, nets, gainA, gainB);
    } while (totalGain > 0 && --pass > 0);


    // cutsize
    cutSize = getCutSize(nets, cells);
    cout << "CutSize: " << cutSize << "\n";

    processOutput(outputFile, dieA, dieB, cutSize);

    auto duration = duration_cast<milliseconds>(high_resolution_clock::now() - start);
    cout << "Execution Time: " << duration.count() / (1000 * 60) << " min, " 
         << (duration.count() % (1000 * 60)) / 1000 << " sec\n";

    return 0;
}