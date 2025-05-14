#include "Net.hpp"

Net::Net() {
    pins = std::vector<Pin*>();
    hardblocks = std::vector<Hardblock*>();
}

void Net::add_pin(Pin* pin) {
    pins.push_back(pin);
}

void Net::add_hardblock(Hardblock* hardblock) {
    hardblocks.push_back(hardblock);
}

int Net::HPWL() {
    int min_x = INT_MAX, max_x = INT_MIN, min_y = INT_MAX, max_y = INT_MIN;
    // std::cout << "Net HPWL:" << std::endl;
    // std::cout << "Pins:" << pins.size() << std::endl;
    // std::cout << "Hardblocks:" << hardblocks.size() << std::endl;
    for(auto pin : pins) {
        // std::cout << "Pin: " << pin->name << std::endl;
        Coord c = pin->coord;
        min_x = std::min(min_x, c.x);
        max_x = std::max(max_x, c.x);
        min_y = std::min(min_y, c.y);
        max_y = std::max(max_y, c.y);
    }
    for(auto hardblock : hardblocks) {
        // std::cout << "Hardblock: " << hardblock->name << std::endl;
        Coord c = hardblock->get_center();
        min_x = std::min(min_x, c.x);
        max_x = std::max(max_x, c.x);
        min_y = std::min(min_y, c.y);
        max_y = std::max(max_y, c.y);
    }
    return (max_x - min_x) + (max_y - min_y);
}