#pragma once
#include <bits/stdc++.h>
#include "Coord.hpp"

struct Hardblock {
    std::string name;
    Coord coord;
    int width, height;
    bool rotated;
    Hardblock():name(),coord(),width(0),height(0),rotated(false){}
    Hardblock(const std::string &n, const Coord &c, int w, int h)
        : name(n), coord(c), width(w), height(h), rotated(false) {}
};