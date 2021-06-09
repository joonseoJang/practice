#include <iostream>
#include "level_set_extractor.h"

int main(int, char**) {
    std::vector<std::vector<double>> polygon = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    levelSetExtractor lse;
    lse.initialize(polygon, 0.1, 0.1, 15, 0.5);
}
