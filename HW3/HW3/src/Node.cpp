#include "Node.hpp"

Node::Node() {
    type = "default";
    index = -1;
    width = -1, height = -1;
    left_from = -1, left_at = -1;
    right_from = -1, right_at = -1;
    coord = Coord(0, 0);
}

Node::Node(std::string type, int index, int width, int height, int left_from, int left_at, int right_from, int right_at, Coord coord) {
    this->type = type;
    this->index = index;
    this->width = width, this->height = height;
    this->left_from = left_from, this->left_at = left_at;
    this->right_from = right_from, this->right_at = right_at;
    this->coord = coord;
}
