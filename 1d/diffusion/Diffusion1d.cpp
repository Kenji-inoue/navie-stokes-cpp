#include <iostream>
#include <string>
#include <stdexcept>
#include "Diffusion1d.h"

Diffusion1d::Diffusion1d(int meshX, double constA, double deltaX, double deltaT) 
    : MESH_X(meshX), CONST_A(constA), DELTA_X(deltaX), DELTA_T(deltaT)
{
    f.resize(MESH_X);
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
        if (4 <= i && i <= 5) {
            f[i] = 1.0;
        }
        else {
            f[i] = 0.0;
        }
        f_next[i] = 0.0;
    }
}

void Diffusion1d::simulate() {
    const int TOTAL_ITERATION = 10;
    for (int time = 0; time < TOTAL_ITERATION; time++) {
        displayProcess(time);
        calculate();
        updateField();
    }
}

void Diffusion1d::displayProcess(int time) {
    const int DISPLAY_INTERVAL = 1;
    if (time % DISPLAY_INTERVAL == 0) {
        printf("t:%d", time);
        for (int i = 0; i < MESH_X; i++) {
            printf("%6.3f", f[i]);
        }
        printf("\n");
    }
}

void Diffusion1d::calculate() {
    const auto coeff = (CONST_A * DELTA_T / DELTA_X / DELTA_X);
    for (int i = 1; i <= MESH_X - 2; i++) {
        f_next[i] = f[i] + coeff * (f[i + 1] - 2 * f[i] + f[i - 1]);
    }
}

void Diffusion1d::updateField() {
    f = f_next;
}
