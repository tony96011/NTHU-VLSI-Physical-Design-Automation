#include "Floorplan.hpp"
#include "Coord.hpp"
#include <bits/stdc++.h>

void Floorplan::setRandomSeed(uint64_t seed) {
    rng.seed(seed);
    srand((unsigned)seed);         // 如果还想保留 rand() 的话
}

Floorplan::Floorplan(const std::string &inFile, const std::string &outFile) {
    // int num_sa = 1;
    // Set random seed using fixed value
    std::cout << "inFile: " << inFile << "\n";
    if(inFile == "../testcase/public1.txt")
        setRandomSeed(1);  // Use a fixed seed value
    else if(inFile == "../testcase/public2.txt"){
        setRandomSeed(2);  // Use a fixed seed value
    }
    else if(inFile == "../testcase/public3.txt")
        setRandomSeed(3);  // Use a fixed seed value
    
    parseInput(inFile);
    buildSymmetryIslands();

    auto [initial_pre, initial_in] = init_sol();
    BStarTree<int64_t> global;
    global.buildTree(initial_pre, initial_in);
    global.setPosition();
    auto [W,H] = global.getWidthHeight(global.root);
    std::cout << "\nInitial area: " << W << "x" << H  << " = " << global.getArea() << "\n";



    // Perform SA twice using a loop
    std::vector<std::string> empty = {};
    auto curr_pre = initial_pre;
    auto curr_in = initial_in;
    
    for(int i = 0; i < 1; i++) {
        auto [new_pre, new_in] = simulatedAnnealing(curr_pre, curr_in, empty);
        
        global.buildTree(new_pre, new_in);
        global.setPosition();
        W = global.getWidthHeight(global.root).first;
        H = global.getWidthHeight(global.root).second;
        std::cout << "\nSA iteration " << i + 1 << " area: " << W << "x" << H << " = " << W*H << "\n";
        
        curr_pre = new_pre;
        curr_in = new_in;
    }

    syncCoordinates(curr_pre);
    int64_t max_W, max_H;
    updateBoundary(max_W, max_H);
    std::cout << "\nFinal area: " << max_W << "x" << max_H  << " = " << max_H*max_W << "\n";
    writeOutput(outFile, max_W, max_H);
    checkOverlap();
}


void Floorplan::updateBoundary(int64_t &totalW, int64_t &totalH) {
    totalW = 0;
    totalH = 0;
    for (auto &pr : HB) {
        auto &b = pr.second;
        int64_t w = b.rotated ? b.height : b.width;
        int64_t h = b.rotated ? b.width  : b.height;
        totalW = std::max(totalW, b.coord.x + w);
        totalH = std::max(totalH, b.coord.y + h);
    }
}


void Floorplan::greedyVerticalPack() {
    // 1) 收集所有 symmetry‐block id
    std::unordered_set<std::string> symIds;
    for (auto &kv : SG) {
        for (auto &p : kv.second.getSymPairs()) {
            symIds.insert(p.first);
            symIds.insert(p.second);
        }
        for (auto &s : kv.second.getSymSelfs())
            symIds.insert(s);
    }

    // 2) 收集所有要做 packing 的非对称块并按 x 排序（从左往右）
    std::vector<Hardblock*> movers;
    movers.reserve(HB.size());
    for (auto &pr : HB) {
        if (!symIds.count(pr.first))
            movers.push_back(&pr.second);
    }
    std::sort(movers.begin(), movers.end(),
              [](auto *a, auto *b){ return a->coord.x < b->coord.x; });

    // 3) 对每个非对称块向下滑动
    for (auto *b : movers) {
        int64_t x0 = b->coord.x;
        int64_t x1 = x0 + (b->rotated ? b->height : b->width);
        int64_t y_old = b->coord.y;

        // 寻找所有与 b 在 x 轴上重叠且位于 b 之下的 "挡路"块的顶边缘
        int64_t bestY = 0;
        for (auto &pr : HB) {
            auto *o = &pr.second;
            if (o == b) continue;
            // x 重叠判定
            int64_t ox0 = o->coord.x;
            int64_t ox1 = ox0 + (o->rotated ? o->height : o->width);
            if (ox1 <= x0 || ox0 >= x1) continue;
            // o 的顶边
            int64_t oy1 = o->coord.y + (o->rotated ? o->width : o->height);
            // 只考虑 o 在 b 当前 y 之下，才能当作阻挡
            if (oy1 <= y_old) {
                bestY = std::max(bestY, oy1);
            }
        }

        // 把 b 滑到 bestY
        b->coord.y = bestY;
    }
}


void Floorplan::greedyHorizontalPack() {
    // 1) 收集所有 symmetry‐block id
    std::unordered_set<std::string> symIds;
    for (auto &kv : SG) {
        for (auto &p : kv.second.getSymPairs()) {
            symIds.insert(p.first);
            symIds.insert(p.second);
        }
        for (auto &s : kv.second.getSymSelfs())
            symIds.insert(s);
    }

    // 2) 收集所有要做 packing 的非对称块并按 y 排序（从下往上）
    std::vector<Hardblock*> movers;
    movers.reserve(HB.size());
    for (auto &pr : HB) {
        if (!symIds.count(pr.first))
            movers.push_back(&pr.second);
    }
    std::sort(movers.begin(), movers.end(),
              [](auto *a, auto *b){ return a->coord.y < b->coord.y; });

    // 3) 对每个非对称块向左滑动
    for (auto *b : movers) {
        int64_t y0 = b->coord.y;
        int64_t y1 = y0 + (b->rotated ? b->width  : b->height);
        int64_t x_old = b->coord.x;

        // 寻找所有与 b 在 y 轴上重叠"挡路"块的右边缘
        // （注意：这里不跳过对称块，让它们也当障碍物）
        int64_t bestX = 0;
        for (auto &pr : HB) {
            auto *o = &pr.second;
            if (o == b) continue;
            // y 方向必须有交叠
            int64_t oy0 = o->coord.y;
            int64_t oy1 = oy0 + (o->rotated ? o->width : o->height);
            if (oy1 <= y0 || oy0 >= y1) continue;
            // o 在 b 左侧时才能挡道
            int64_t ox1 = o->coord.x + (o->rotated ? o->height : o->width);
            if (ox1 <= x_old) {
                bestX = std::max(bestX, ox1);
            }
        }

        // 把 b 滑到 bestX
        b->coord.x = bestX;
    }
}


void Floorplan::checkOverlap() {
    bool hasOverlap = false;
    for (auto &pr : HB) {
        auto &hb = pr.second;
        for (auto &pr2 : HB) {
            auto &hb2 = pr2.second;
            if (hb.name != hb2.name) {
                // Get the actual dimensions considering rotation
                int64_t hb_width = hb.rotated ? hb.height : hb.width;
                int64_t hb_height = hb.rotated ? hb.width : hb.height;
                int64_t hb2_width = hb2.rotated ? hb2.height : hb2.width;
                int64_t hb2_height = hb2.rotated ? hb2.width : hb2.height;

                // Check if two rectangles overlap
                bool xOverlap = (hb.coord.x < hb2.coord.x + hb2_width) && 
                              (hb.coord.x + hb_width > hb2.coord.x);
                bool yOverlap = (hb.coord.y < hb2.coord.y + hb2_height) && 
                              (hb.coord.y + hb_height > hb2.coord.y);
                
                if (xOverlap && yOverlap) {
                    hasOverlap = true;
                    std::cout << "Overlap detected between " << hb.name << " and " << hb2.name << "\n";
                    std::cout << hb.name << ": (" << hb.coord.x << "," << hb.coord.y 
                             << ") " << hb_width << "x" << hb_height 
                             << (hb.rotated ? " (rotated)" : "") << "\n";
                    std::cout << hb2.name << ": (" << hb2.coord.x << "," << hb2.coord.y 
                             << ") " << hb2_width << "x" << hb2_height
                             << (hb2.rotated ? " (rotated)" : "") << "\n";
                }
            }
        }
    }
    if (!hasOverlap) {
        std::cout << "No overlaps detected.\n";
    }
}

void Floorplan::syncCoordinates(std::vector<Node<int64_t>*>& pre) {
    for(auto &node: pre) {
        // std::cout << "node: " << node->id << "\n";
        if(node->isHierNode()) {
            auto *hnode = static_cast<HierNode<int64_t>*>(node);
            auto *local_node = hnode->local_tree.root;
            // std::cout << "hier_node: " << hnode->id << "\n";
            // std::cout << "hnode->x: " << hnode->x << "\n";
            // std::cout << "hnode->y: " << hnode->y << "\n";
            // std::cout << "hnode->width: " << hnode->width << "\n";
            // std::cout << "hnode->height: " << hnode->height << "\n";
            
            auto nodes = preorderTraversal(local_node, hnode->size);

            for(auto &n: nodes) {
                if(HB.find(n->id) != HB.end()) {
                    // std::cout << "n->id: " << n->id << "\n";
                    // std::cout << "n->x: " << n->x << "\n";
                    // std::cout << "n->y: " << n->y << "\n";
                    // std::cout << "n->width: " << n->width << "\n";
                    // std::cout << "n->height: " << n->height << "\n";

                    // check if the block has a symmetric block
                    if(hnode->sym_map.find(n->id) != hnode->sym_map.end()) {
                        HB[n->id].coord.x = n->x + hnode->x + hnode->width/2;
                        HB[n->id].coord.y = n->y + hnode->y;
                        if(n->width != HB[n->id].width || n->height != HB[n->id].height) {
                            HB[n->id].rotated = true;
                        }
                        auto *sym_node = hnode->sym_map[n->id];
                        HB[sym_node->id].coord.x = hnode->width/2 - (n->x + n->width) + hnode->x;
                        HB[sym_node->id].coord.y = n->y + hnode->y;
                        HB[sym_node->id].rotated = HB[n->id].rotated;
                    }
                    else{
                        HB[n->id].coord.x = hnode->x + hnode->width/2 - n->width;
                        HB[n->id].coord.y = n->y + hnode->y;
                        if(n->width != HB[n->id].width && n->height != HB[n->id].height) {
                            HB[n->id].rotated = true;
                        }
                    }
                }
            }
        } else {
            if(HB.find(node->id) != HB.end()) {
                HB[node->id].coord.x = node->x;
                HB[node->id].coord.y = node->y;
                // Check if the block is rotated
                if(node->width != HB[node->id].width || node->height != HB[node->id].height) {
                    HB[node->id].rotated = true;
                }
            }
        }
    }
}

void Floorplan::parseInput(const std::string &f) {
    std::ifstream fin(f);
    if (!fin) { std::cerr<<"Cannot open "<<f<<"\n"; std::exit(1); }
    std::string line, key;
    while (std::getline(fin,line)) {
        if (line.empty()||line.rfind("//",0)==0) continue;
        std::istringstream is(line);
        is >> key;
        if (key=="NumHardBlocks") {
            is >> NH;
            for (int i=0; i<NH; ) {
                std::getline(fin,line);
                if (line.empty()||line.rfind("//",0)==0) continue;
                std::istringstream bs(line);
                std::string t,n; int w,h;
                bs >> t >> n >> w >> h;
                HB[n] = Hardblock(n, Coord(0,0), w, h);
                // std::cout << n << " " << w << " " << h << "\n";
                auto *node = new Node<int64_t>(n, w, h);
                nodes[n] = node;
                hier.push_back(n);
                ++i;
            }
        }
        else if (key=="NumSymGroups") {
            is >> NS;
            for (int i=0; i<NS; ) {
                std::getline(fin,line);
                if (line.empty()||line.rfind("//",0)==0) continue;
                std::istringstream gs(line);
                std::string t,g; int cnt;
                gs >> t >> g >> cnt;
                SymGroup grp;
                hier.push_back(g);
                for (int j=0; j<cnt; ) {
                    std::getline(fin,line);
                    if (line.empty()||line.rfind("//",0)==0) continue;
                    std::istringstream ss(line);
                    std::string et; ss>>et;
                    if (et=="SymPair") {
                        std::string a,b; ss>>a>>b;
                        grp.addSymPair(a,b);
                        hier.erase(std::remove(hier.begin(), hier.end(), a), hier.end());
                        hier.erase(std::remove(hier.begin(), hier.end(), b), hier.end());
                    } else {
                        std::string a; ss>>a;
                        grp.addSymSelf(a);
                        hier.erase(std::remove(hier.begin(), hier.end(), a), hier.end());
                    }
                    ++j;
                }
                SG[g] = grp;
                auto *hnode = new HierNode<int64_t>(g, grp.getSymPairs(), grp.getSymSelfs());
                nodes[g] = hnode;
                ++i;
            }
        }
    }
}


void Floorplan::buildSymmetryIslands() {
    for(auto &node: nodes){
        auto *n = node.second;
        if(n->isHierNode()){
            auto *hnode = static_cast<HierNode<int64_t>*>(n);
            BStarTree<int64_t> localTree = packIsland(hnode);
            hnode->local_tree = std::move(localTree);
        }
    }

    for(auto &name: hier){
        nodes_hier[name] = nodes[name];
    }
}


BStarTree<int64_t> Floorplan::packIsland(HierNode<int64_t> *hnode) {
    auto &pairs = hnode->pairs;
    auto &selfs = hnode->selfs;

    // Clear the sym_map before populating it
    hnode->sym_map.clear();

    // First handle selfs
    for (int i=0; i<(int)selfs.size()-1; ++i) {
        if (nodes.find(selfs[i]) == nodes.end() || nodes.find(selfs[i+1]) == nodes.end()) {
            std::cout << "Warning: self node not found in nodes map\n";
            continue;
        }
        auto *n = nodes[selfs[i]];
        auto *m = nodes[selfs[i+1]];
        if (!n || !m) {
            std::cout << "Warning: null node encountered in selfs\n";
            continue;
        }
        n->width = n->width/2;
        n->rchild = m;
        m->parent = n;
    }

    if(selfs.size() == 1){
        auto *n = nodes[selfs[0]];
        n->width = n->width/2;
        n->rchild = nullptr;
        n->parent = nullptr;
    }

    // Then handle pairs
    for(int i=0; i<(int)pairs.size()-1; ++i) {
        if (nodes.find(pairs[i].first) == nodes.end() || nodes.find(pairs[i+1].first) == nodes.end()) {
            std::cout << "Warning: pair node not found in nodes map\n";
            continue;
        }
        auto *n = nodes[pairs[i].first];
        auto *m = nodes[pairs[i+1].first];
        if (!n || !m) {
            std::cout << "Warning: null node encountered in pairs\n";
            continue;
        }
        n->rchild = m;
        m->parent = n;
    }

    // Connect selfs and pairs if both exist
    if(selfs.size() > 0 && pairs.size() > 0) {
        if (nodes.find(selfs[0]) == nodes.end() || nodes.find(pairs[0].first) == nodes.end()) {
            std::cout << "Warning: node not found when connecting selfs and pairs\n";
        } else {
            auto *n = nodes[selfs[0]];
            auto *m = nodes[pairs[0].first];
            if (n && m) {
                n->lchild = m;
                m->parent = n;
            }
        }
    }

    // Populate sym_map for pairs
    for(auto &pair: pairs) {
        if (nodes.find(pair.first) == nodes.end() || nodes.find(pair.second) == nodes.end()) {
            std::cout << "Warning: symmetric pair node not found in nodes map\n";
            continue;
        }
        auto *first = nodes[pair.first];
        auto *second = nodes[pair.second];
        if (first && second) {
            hnode->sym_map[pair.first] = second;  // Map first to second
            hnode->sym_map[pair.second] = first;  // Also map second to first for completeness
        }
    }

    // Calculate size and get root
    int K = (int)selfs.size() + (int)pairs.size();
    hnode->size = K;
    Node<int64_t> *root = nullptr;
    
    if((int)selfs.size() > 0) {
        if (nodes.find(selfs[0]) != nodes.end()) {
            root = nodes[selfs[0]];
            // std::cout << "Root set to first self node: " << selfs[0] << "\n";
        } else {
            std::cout << "Warning: First self node not found in nodes map\n";
        }
    } else if((int)pairs.size() > 0) {
        if (nodes.find(pairs[0].first) != nodes.end()) {
            root = nodes[pairs[0].first];
            // std::cout << "Root set to first pair node: " << pairs[0].first << "\n";
        } else {
            std::cout << "Warning: First pair node not found in nodes map\n";
        }
    }

    if (!root) {
        std::cout << "Error: Could not find root node. Selfs size: " << selfs.size() 
                  << ", Pairs size: " << pairs.size() << "\n";
        return BStarTree<int64_t>();
    }
    
    std::vector<Node<int64_t>*> pre = preorderTraversal(root, K);
    std::vector<Node<int64_t>*> in  = inorderTraversal(root, K);

    if (pre.empty() || in.empty()) {
        std::cout << "Error: Empty traversal results\n";
        return BStarTree<int64_t>();
    }

    // use simulated annealing to find the tree with minimum area
    auto [best_pre, best_in] = simulatedAnnealing(pre, in, selfs);

    BStarTree<int64_t> tree;
    tree.buildTree(best_pre,best_in);
    // tree.buildTree(pre, in);
    tree.setPosition();
    auto [W, H] = tree.getWidthHeight(tree.root);
    hnode->width = 2*W; hnode->height = H;

    return tree;
}

std::pair<std::vector<Node<int64_t>*>,std::vector<Node<int64_t>*>> Floorplan::init_sol() {
    if (nodes_hier.empty()) return {{},{}};

    for (auto &kv : nodes_hier) {
        kv.second->lchild = kv.second->rchild = nullptr;
        kv.second->parent = nullptr;
    }

    auto parent_it = nodes_hier.begin();
    auto child_it  = std::next(parent_it);

    Node<int64_t>* root = parent_it->second;

    while (parent_it != nodes_hier.end() && child_it != nodes_hier.end()) {
        Node<int64_t>* u = parent_it->second;
        u->lchild = child_it->second;
        child_it->second->parent = u;
        ++child_it;
        if (child_it != nodes_hier.end()) {
            u->rchild = child_it->second;
            child_it->second->parent = u;
            ++child_it;
        }
        ++parent_it;
    }

    size_t M = nodes_hier.size();
    auto pre = preorderTraversal(root, (int)M);
    auto in  = inorderTraversal (root, (int)M);

    return { std::move(pre), std::move(in) };
}

void Floorplan::move_node(std::vector<Node<int64_t>*>& pre, std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs) {
    if (pre.size() <= 1) return;  // Need at least 2 nodes to move
    
    // Find a random non-root node to move
    int idx = 1 + rand() % (pre.size() - 1);  // Skip root (index 0)
    Node<int64_t>* node = pre[idx];
    
    // Skip if node is a HierNode
    if (node->isHierNode()) return;
    
    // Find a random target position
    int target_idx = rand() % pre.size();
    while(target_idx == idx) {
        target_idx = rand() % pre.size();
    }
    Node<int64_t>* target = pre[target_idx];
    
    // Skip if target is a HierNode
    if (target->isHierNode()) return;
    
    // Skip if target is a descendant of node
    Node<int64_t>* current = target;
    while (current && current != node) {
        current = current->parent;
    }
    if (current == node) return;  // Target is a descendant of node
    
    // Remove node from its current position
    if (node->parent) {
        if (node->parent->lchild == node) {
            node->parent->lchild = nullptr;
        } else {
            node->parent->rchild = nullptr;
        }
    }
    
    // Save node's children
    Node<int64_t>* left_child = node->lchild;
    Node<int64_t>* right_child = node->rchild;
    
    // Clear node's children and parent
    node->lchild = nullptr;
    node->rchild = nullptr;
    node->parent = nullptr;
    
    // Insert node as a child of target
    if (rand() % 2) {  // Insert as left child
        if (target->lchild) {
            // If target already has a left child, make it the right child of our node
            node->rchild = target->lchild;
            target->lchild->parent = node;
        }
        target->lchild = node;
        node->parent = target;
    } else {  // Insert as right child
        if (target->rchild) {
            // If target already has a right child, make it the left child of our node
            node->lchild = target->rchild;
            target->rchild->parent = node;
        }
        target->rchild = node;
        node->parent = target;
    }
    
    // Reattach the original children
    if (left_child) {
        if (!node->lchild) {
            node->lchild = left_child;
            left_child->parent = node;
        } else {
            // Find the rightmost node in the left subtree
            Node<int64_t>* current = node->lchild;
            while (current->rchild) {
                current = current->rchild;
            }
            current->rchild = left_child;
            left_child->parent = current;
        }
    }
    
    if (right_child) {
        if (!node->rchild) {
            node->rchild = right_child;
            right_child->parent = node;
        } else {
            // Find the leftmost node in the right subtree
            Node<int64_t>* current = node->rchild;
            while (current->lchild) {
                current = current->lchild;
            }
            current->lchild = right_child;
            right_child->parent = current;
        }
    }
    
    // Update traversals
    int n = pre.size();
    pre = preorderTraversal(pre[0], n);
    in = inorderTraversal(pre[0], n);
    
    // Verify the number of nodes
    if (pre.size() != n) {
        std::cout << "Warning: Node count changed from " << n << " to " << pre.size() << " in move_node\n";
        // Restore the original traversals
        pre = preorderTraversal(pre[0], n);
        in = inorderTraversal(pre[0], n);
    }
}

void Floorplan::perturb(std::vector<Node<int64_t>*>& pre,
    std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs) {
    double r = static_cast<double>(rand()) / RAND_MAX;
    if (r < 0.2) {
        rotate_node(pre, selfs);
    } else  if (r < 0.6) {
        swap_node(pre, in, selfs);
    } else {
        move_node(pre, in, selfs);
    }
}

void Floorplan::swap_node(std::vector<Node<int64_t>*>& pre, std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs) {
    if (pre.empty()) return;
    
    // Find two random nodes
    int idx1 = rand() % pre.size();
    Node<int64_t>* node1 = pre[idx1];

    int idx2 = rand() % pre.size();
    while(idx1 == idx2) {
        idx2 = rand() % pre.size();
    }
    Node<int64_t>* node2 = pre[idx2];

    //if selected node is a self symmetry node, return
    if(std::find(selfs.begin(), selfs.end(), node1->id) != selfs.end() || std::find(selfs.begin(), selfs.end(), node2->id) != selfs.end()) return;
    
    // Find their positions in inorder traversal
    auto it1 = std::find(in.begin(), in.end(), node1);
    auto it2 = std::find(in.begin(), in.end(), node2);
    if (it1 == in.end() || it2 == in.end()) return;

    // Swap in preorder
    std::swap(pre[idx1], pre[idx2]);

    // Swap in inorder
    std::swap(*it1, *it2);
}


void Floorplan::rotate_node(std::vector<Node<int64_t>*>& pre, std::vector<std::string>& selfs) {
    int idx = rand() % pre.size();
    auto *node = pre[idx];
    if(!node->isHierNode()) {
        if(std::find(selfs.begin(), selfs.end(), node->id) != selfs.end()){
            node->width = node->width*2;
            std::swap(node->width, node->height);
            node->width = node->width/2;
        }
        else std::swap(node->width, node->height);
    }
}

std::pair<std::vector<Node<int64_t>*>,std::vector<Node<int64_t>*>> Floorplan::clone_traversals(const std::vector<Node<int64_t>*>& pre, const std::vector<Node<int64_t>*>& in) 
{
    std::vector<Node<int64_t>*> newPre;
    std::vector<Node<int64_t>*> newIn;
    std::unordered_map<std::string, Node<int64_t>*> id_to_node;

    // First pass: create all nodes
    for(auto &node : pre) {
        if (!node) continue;
        
        Node<int64_t>* newNode;
        if (node->isHierNode()) {
            auto *old_hnode = static_cast<HierNode<int64_t>*>(node);
            auto *new_hnode = new HierNode<int64_t>(old_hnode->id, old_hnode->pairs, old_hnode->selfs);
            new_hnode->width = old_hnode->width;
            new_hnode->height = old_hnode->height;
            new_hnode->size = old_hnode->size;
            new_hnode->sym_map = old_hnode->sym_map;
            
            // Clone local_tree
            if (old_hnode->local_tree.root) {
                auto local_pre = preorderTraversal(old_hnode->local_tree.root, old_hnode->size);
                auto local_in = inorderTraversal(old_hnode->local_tree.root, old_hnode->size);
                new_hnode->local_tree.buildTree(local_pre, local_in);
                new_hnode->local_tree.setPosition();
            }
            
            newNode = new_hnode;
        } else {
            newNode = new Node<int64_t>(node->id, node->width, node->height);
        }
        newPre.push_back(newNode);
        id_to_node[node->id] = newNode;
    }

    // Second pass: create inorder traversal
    for(auto &node : in) {
        if (!node) continue;
        auto it = id_to_node.find(node->id);
        if(it != id_to_node.end()) {
            newIn.push_back(it->second);
        } else {
            std::cerr << "Error: Node " << node->id << " not found in newPre\n";
        }
    }
    
    return {newPre, newIn};
}


std::pair<std::vector<Node<int64_t>*>, std::vector<Node<int64_t>*>>
Floorplan::simulatedAnnealing(
    const std::vector<Node<int64_t>*>& pre,
    const std::vector<Node<int64_t>*>& in, std::vector<std::string>& selfs
) {
    // SA parameters
    double T      = 1000.0;
    const double T_MIN   = 0.01;
    const double T_DECAY = 0.90;
    const int K       = 10;
    const int N       = static_cast<int>(pre.size()) * K;
    const int DOUBLE_N = N * 2;

    // Deep clone original traversals
    auto [currPre, currIn] = clone_traversals(pre, in);
    auto bestPre = currPre;
    auto bestIn  = currIn;

    // Build initial best tree to compute area
    BStarTree<int64_t> tree;
    tree.buildTree(currPre, currIn);
    tree.setPosition();
    auto [w0, h0] = tree.getWidthHeight(tree.root);
    double bestArea    = static_cast<double>(w0) * h0;
    double currentArea = bestArea;

    // Main SA loop
    while (T > T_MIN) {
        for (int i = 0; i < DOUBLE_N; ++i) {
            // Deep clone current solution
            auto [candPre, candIn] = clone_traversals(currPre, currIn);
            
            BStarTree<int64_t> candidate;
            candidate.buildTree(candPre, candIn);

            perturb(candPre, candIn, selfs);

            candidate.buildTree(candPre, candIn);
            candidate.setPosition();
            auto [w, h] = candidate.getWidthHeight(candidate.root);
            double candArea = static_cast<double>(w) * h;

            double delta = candArea - currentArea;
            if (delta < 0 || std::exp(-delta / T) > static_cast<double>(rand()) / RAND_MAX) {
                // Accept candidate
                currPre = std::move(candPre);
                currIn = std::move(candIn);
                currentArea = candArea;

                // Update best if improved
                if (currentArea < bestArea) {
                    // Free old best solution
                    for (auto* node : bestPre) {
                        if (node) delete node;
                    }
                    
                    bestArea = currentArea;
                    bestPre = currPre;
                    bestIn = currIn;
                }
            } else {
                // Free rejected candidate
                for (auto* node : candPre) {
                    if (node) delete node;
                }
               
            }
        }
        // Cool down
        T *= T_DECAY;
    }

    // Store the final best solution back to nodes map
    for (auto* node : bestPre) {
        if (nodes.find(node->id) != nodes.end()) {
            delete nodes[node->id];  // Delete old node
        }
        nodes[node->id] = node;  // Store new node
    }

    return {bestPre, bestIn};
}

void Floorplan::writeOutput(const std::string &f, int64_t W, int64_t H) {
    std::ofstream fo(f);
    auto area = W * H;
    fo << "Area " << area << "\n";
    fo<<"NumHardBlocks "<<NH<<"\n";
    for (auto &pr : HB) {
        auto &hb = pr.second;
        fo<<hb.name<<" "<<hb.coord.x<<" "<<hb.coord.y<<" "<<hb.rotated<<"\n";
    }
}