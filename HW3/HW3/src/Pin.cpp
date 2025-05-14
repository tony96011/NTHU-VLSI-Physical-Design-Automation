#include "Pin.hpp"

Pin::Pin() {
    name = "default";
    coord = Coord();
}

Pin::Pin(std::string n, Coord coord) {
    name = n;
    this->coord = coord;
}