#include "Diffusion1d.h"
#include <stdexcept>
#include <iostream>

int main() {
    try {
        Diffusion1d diffusion(10, 1e-5, 0.001, 0.01);
        diffusion.simulate();
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}