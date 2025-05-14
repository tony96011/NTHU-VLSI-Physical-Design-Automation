#include "Hardblock.hpp"

Hardblock::Hardblock() {
    name = "default";
    coord = Coord();
    width = -1;
    height = -1;
    rotated = -1;
}

Hardblock::Hardblock(std::string name, Coord coord, int width, int height) {
    this->name = name;
    this->coord = coord;
    this->width = width;
    this->height = height;
    rotated = 0;
}

Coord Hardblock::get_center() {
    if(rotated)
        return Coord(coord.x + height / 2, coord.y + width / 2);
    return Coord(coord.x + width / 2, coord.y + height / 2);
}
