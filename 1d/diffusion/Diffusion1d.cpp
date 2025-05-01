#include <iostream>
#include <string>
#include <stdexcept>
#include "Diffusion1d.h"

Diffusion1d::Diffusion1d(int meshX, double constA, double deltaX, double deltaT) 
    : MESH_X(meshX), CONST_A(constA), DELTA_X(deltaX), DELTA_T(deltaT)
{
    f_next.resize(MESH_X);
    validateTime();
    initializeField();
}

void Diffusion1d::validateTime() {
    const auto DELTA_T_MAX = 0.2 * DELTA_X * DELTA_X / CONST_A;
    if(DELTA_T > DELTA_T_MAX) {
        throw std::runtime_error("DELTA_T is too large. Must be <= " + std::to_string(DELTA_T_MAX));
    }
}

void Diffusion1d::initializeField() {
    for (int i = 0; i < MESH_X; i++) {
        f_next[i] = 0.0;
    }
}

void Diffusion1d::displayProcess(int time, int interval, const Field1d& f) {
    if (time % interval != 0) {
        return;
    }

    printf("t:%d", time);
    for (int i = 0; i < MESH_X; i++) {
        printf("%6.3f", f[i]);
    }
    printf("\n");
}

Field1d Diffusion1d::calculate(const Field1d& f) {
    const auto coeff = (CONST_A * DELTA_T / DELTA_X / DELTA_X);
    for (int i = 1; i <= MESH_X - 2; i++) {
        f_next[i] = f[i] + coeff * (f[i + 1] - 2 * f[i] + f[i - 1]);
    }
    return f_next;
}
