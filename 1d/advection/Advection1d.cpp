#include <iostream>
#include <string>
#include <stdexcept>
#include "Advection1d.h"

Advection1d::Advection1d(int meshX, double constU, double deltaX, double deltaT) 
    : MESH_X(meshX), CONST_U(constU), DELTA_X(deltaX), DELTA_T(deltaT)
{
    validateTime();
}

void Advection1d::validateTime() {
    const auto DELTA_T_MAX = 0.2 * DELTA_X / CONST_U;
    if(DELTA_T > DELTA_T_MAX) {
        throw std::runtime_error("DELTA_T is too large. Must be <= " + std::to_string(DELTA_T_MAX));
    }
}

void Advection1d::displayProcess(int time, int interval, const Field1d& f) {
    if (time % interval != 0) {
        return;
    }

    printf("t:%d", time);
    for (int i = 0; i < MESH_X; i++) {
        printf("%6.3f", f[i]);
    }
    printf("\n");
}

Field1d Advection1d::calculate(const Field1d& f) {
    Field1d f_next(MESH_X);
    for (int i = 1; i <= MESH_X - 2; i++) {
        if (CONST_U >= 0) {
            f_next[i] = f[i] + (CONST_U * DELTA_T / DELTA_X) * (f[i - 1] - f[i]);
        }
        else {
            f_next[i] = f[i] + (CONST_U * DELTA_T / DELTA_X) * (f[i] - f[i+1]);
        }
    }
    periodicBoundaryCondition(f_next);
    return f_next;
}

void Advection1d::periodicBoundaryCondition(Field1d& f) {
    f[0] = f[MESH_X - 2];
    f[MESH_X - 1] = f[1];
}
