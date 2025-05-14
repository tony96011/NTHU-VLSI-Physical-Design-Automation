# B*-Tree Usage
1. Create a Node vector as the initial solution.
2. Pass the Node vector into the function `buildTree` to build a B*-Tree.
3. Pass the root node returned from the function `buildTree` into the function `setPosition`.
4. After setting the position of all the nodes, call the function `getArea` to get the area of the placement result.

Example:
```cpp
std::vector<Node<int64_t> *> initSolution = yourOwnSolution();
BStarTree<int64_t> bStarTree;
Node<int64_t> *root = bStarTree.buildTree(initSolution);
bStarTree.setPosition(root);
int64_t area = bStarTree.getArea(root);
```
