// BStarTree.hpp
#include "SegmentTree.hpp"
#include <bits/stdc++.h>

/**
 * @brief Node in the B*-tree for floorplanning.
 */
template <typename T>
struct Node
{
    using ptr = std::unique_ptr<Node>;
    std::string id;
    T x, y;
    T width, height;
    Node *lchild, *rchild, *parent;

    Node() : x(0), y(0), width(0), height(0), lchild(nullptr), rchild(nullptr), parent(nullptr) {}
    Node(const std::string &id_) : x(0), y(0), width(0), height(0), lchild(nullptr), rchild(nullptr), parent(nullptr), id(id_) {}
    Node(const std::string &id_, T width_, T height_) : x(0), y(0), width(width_), height(height_), lchild(nullptr), rchild(nullptr), parent(nullptr), id(id_) {}
    
    void setPosition(T x_, T y_)
    {
        x = x_;
        y = y_;
    }

    void setShape(T width_, T height_)
    {
        width = width_;
        height = height_;
    }

    bool isHierNode()
    {
        if(id.find("sg") == 0) return true;
        else return false;
    }
};


inline std::vector<Node<int64_t>*> preorderTraversal(Node<int64_t>* root, int size) {
    if (!root) {
        std::cout << "Warning: root is null in preorderTraversal\n";
        return {};
    }
    
    std::vector<Node<int64_t>*> result;
    try {
        result.reserve(size);
    } catch (const std::exception& e) {
        std::cout << "Warning: Failed to reserve space in preorderTraversal: " << e.what() << "\n";
        return {};
    }

    std::function<void(Node<int64_t>*)> dfs = [&](Node<int64_t>* u) {
        if (!u) return;
        result.push_back(u);
        if (u->lchild) dfs(u->lchild);
        if (u->rchild) dfs(u->rchild);
    };

    try {
        dfs(root);
    } catch (const std::exception& e) {
        std::cout << "Warning: Exception in preorderTraversal: " << e.what() << "\n";
        return {};
    }
    
    return result;
}

/// @brief Return the inorder (left–root–right) sequence of the tree.
inline std::vector<Node<int64_t>*> inorderTraversal(Node<int64_t>* root, int size) {
    if (!root) {
        std::cout << "Warning: root is null in inorderTraversal\n";
        return {};
    }
    
    std::vector<Node<int64_t>*> result;
    try {
        result.reserve(size);
    } catch (const std::exception& e) {
        std::cout << "Warning: Failed to reserve space in inorderTraversal: " << e.what() << "\n";
        return {};
    }

    std::function<void(Node<int64_t>*)> dfs = [&](Node<int64_t>* u) {
        if (!u) return;
        if (u->lchild) dfs(u->lchild);
        result.push_back(u);
        if (u->rchild) dfs(u->rchild);
    };

    try {
        dfs(root);
    } catch (const std::exception& e) {
        std::cout << "Warning: Exception in inorderTraversal: " << e.what() << "\n";
        return {};
    }
    
    return result;
}

/**
 * @brief A B*-tree to calculate the coordinates of nodes and the area of placement
 */
template <typename T>
class BStarTree
{
    std::unordered_map<Node<T> *, int64_t> toInorderIdx;

    SegmentTree<T> contourH_up, contourH_low;
    SegmentTree<T> contourV_right, contourV_left;

    Node<T> *buildTree(const std::vector<Node<T> *> &preorder, const std::vector<Node<T> *> &inorder, size_t &i, int64_t l, int64_t r)
    {
        if (l > r || i >= preorder.size())
            return nullptr;

        Node<T> *node = preorder[i++];
        assert(toInorderIdx.count(node) > 0 && "Node not found in inorder map.");
        int64_t idx = toInorderIdx[node];
        
        node->lchild = buildTree(preorder, inorder, i, l, idx - 1);
        if (node->lchild) node->lchild->parent = node;
        
        node->rchild = buildTree(preorder, inorder, i, idx + 1, r);
        if (node->rchild) node->rchild->parent = node;
        
        return node;
    }

    T getTotalWidth(Node<T> *node) const
    {
        if (!node)
            return 0;

        return node->width + getTotalWidth(node->lchild) + getTotalWidth(node->rchild);
    }

    T getTotalHeight(Node<T> *node) const
    {
        if (!node)
            return 0;

        return node->height + getTotalHeight(node->lchild) + getTotalHeight(node->rchild);
    }

    void setPosition(Node<T> *node, T startX)
    {
        if (!node)
            return;

        T endX = startX + node->width;
        T y = contourH_up.query(startX, endX - 1);
        contourH_up.update(startX, endX - 1, y + node->height);
        node->setPosition(startX, y);
        setPosition(node->lchild, endX);
        setPosition(node->rchild, startX);
    }

    // void setPosition_sym(Node<T>* node, T startX, T startY){
    //     if (!node) return;

    //     T endX = startX + node->width;
    //     T y_up  = contourH_up.query(startX, endX - 1);
    //     T y_low = contourH_low.query(startX, endX - 1);
    //     T y = y_up;

    //     contourH_up.update(startX, endX - 1, y + node->height);
    //     contourH_low.update(startX, endX - 1, y); 

    //     T endY = startY + node->height;
    //     T x_right = contourV_right.query(startY, endY - 1);
    //     T x_left  = contourV_left.query(startY, endY - 1);
    //     T x = x_right;
  
    //     contourV_right.update(startY, endY - 1, x + node->width);
    //     contourV_left.update(startY, endY - 1, x);

    //     node->setPosition(x, y);

    //     setPosition_sym(node->lchild, x + node->width, y);
    //     setPosition_sym(node->rchild, x, y + node->height);
    // }

    // void setPosition(Node<T> *node, T startX)
    // {
    //     if (!node)
    //         return;

    //     T endX = startX + node->width;
    //     T y = contourH.query(startX, endX - 1);
    //     contourH.update(startX, endX - 1, y + node->height);
    //     node->setPosition(startX, y);
    //     setPosition(node->lchild, endX);
    //     setPosition(node->rchild, startX);
    // }

public:
    Node<T> *root;

    BStarTree() : root(nullptr) {}

    void buildTree(const std::vector<Node<T> *> &preorder, const std::vector<Node<T> *> &inorder)
    {
        assert(preorder.size() == inorder.size() && "The size of preorder and inorder must be the same.");
        int64_t n = inorder.size();

        toInorderIdx.clear();
        for (int64_t i = 0; i < n; ++i)
            toInorderIdx[inorder[i]] = i;

        size_t i = 0;
        root = buildTree(preorder, inorder, i, 0LL, n - 1);
    }

    void setPosition()
    {
        contourH_up.init(getTotalWidth(root), SegmentTree<T>::Type::MAX);
        setPosition(root, 0);
    }

    // void setPosition_sym()
    // {
    //     contourH_up.init(getTotalWidth(root), SegmentTree<int>::Type::MAX);
    //     contourH_low.init(getTotalWidth(root), SegmentTree<int>::Type::MIN);
    //     contourV_right.init(getTotalHeight(root), SegmentTree<int>::Type::MAX);
    //     contourV_left.init(getTotalHeight(root), SegmentTree<int>::Type::MIN);
    //     setPosition_sym(root, 0, 0);
    // }

    T getArea() const
    {
        auto [width, height] = getWidthHeight(root);
        return width * height;
    }

    std::pair<T, T> getWidthHeight(Node<T> *node) const
    {
        if (!node)
            return {0, 0};

        auto [lMaxWidth, lMaxHeight] = getWidthHeight(node->lchild);
        auto [rMaxWidth, rMaxHeight] = getWidthHeight(node->rchild);
        T maxWidth = std::max(std::max(lMaxWidth, rMaxWidth), node->x + node->width);
        T maxHeight = std::max(std::max(lMaxHeight, rMaxHeight), node->y + node->height);
        return {maxWidth, maxHeight};
    }
};

template<typename T>
struct HierNode : public Node<T> {
    BStarTree<T> local_tree;
    std::vector<std::pair<std::string,std::string>> pairs;
    std::vector<std::string> selfs;
    T size;
    std::unordered_map<std::string, Node<T>*> sym_map;  // Maps surviving node id to its symmetric node pointer
    bool vertical;
    HierNode(const std::string &id_, std::vector<std::pair<std::string,std::string>> pairs, std::vector<std::string> selfs): Node<T>(id_), pairs(pairs), selfs(selfs), size(0), vertical(true) {}
};