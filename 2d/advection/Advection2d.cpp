#include <string>
#include <stdexcept>
#include "Advection2d.h"
#include "FieldUtil.h"

Advection2d::Advection2d(int meshX, int meshY, double constU, double constV, double deltaX, double deltaY, double deltaT) 
    : MESH_X(meshX), MESH_Y(meshY), CONST_U(constU), CONST_V(constV), DELTA_X(deltaX), DELTA_Y(deltaY), DELTA_T(deltaT)
{
    validateTime();
}

void Advection2d::validateTime() {
    double timeCoefficient;
    if (CONST_U == 0 && CONST_V == 0) {
        throw std::runtime_error("Error!! ALL velocity have zero value! Either velocity should be a non-zero value.");
    }
    else if (CONST_U == 0) {
        timeCoefficient = DELTA_Y / CONST_V;
    }
    else if (CONST_V == 0) {
        timeCoefficient = DELTA_X / CONST_U;
    }
    else {
        timeCoefficient = std::min(DELTA_Y / CONST_V, DELTA_X / CONST_U);
    }
    const auto deltaTMax = 0.2 * timeCoefficient;
    if (DELTA_T > deltaTMax) {
        const std::string errorMessage = "Invalid dt error!!:dt must be <= " + std::to_string(deltaTMax) + " from CFL condition";
        throw std::runtime_error(errorMessage);
    }
}

Field2d Advection2d::calculate(const Field2d& f) {
    Field2d f_next;
    FieldUtil::setSize(f_next, MESH_X, MESH_Y);

    for (int j = 1; j <= MESH_Y - 2; j++) {
        for (int i = 1; i <= MESH_X - 2; i++) {
            if (CONST_U >= 0 && CONST_V >= 0) {
                f_next[j][i] = f[j][i] + CONST_U * DELTA_T / DELTA_X * (f[j][i - 1] - f[j][i]) + CONST_V * DELTA_T / DELTA_Y * (f[j-1][i] - f[j][i]);
            }
            else if (CONST_U < 0 && CONST_V >= 0) {
                f_next[j][i] = f[j][i] + CONST_U * DELTA_T / DELTA_X * (f[j][i] - f[j][i+1]) + CONST_V * DELTA_T / DELTA_Y * (f[j-1][i] - f[j][i]);
            }
            else if (CONST_U >= 0 && CONST_V < 0) {
                f_next[j][i] = f[j][i] + CONST_U * DELTA_T / DELTA_X * (f[j][i-1] - f[j][i]) + CONST_V * DELTA_T / DELTA_Y * (f[j][i] - f[j+1][i]);
            }
            else {
                f_next[j][i] = f[j][i] + CONST_U * DELTA_T / DELTA_X * (f[j][i] - f[j][i+1]) + CONST_V * DELTA_T / DELTA_Y * (f[j][i] - f[j+1][i]);
            }
        }
    }
    updateBoundaryCondition(f_next);
    return f_next;
}

void Advection2d::updateBoundaryCondition(Field2d& f) {
    for (int j = 0; j < MESH_Y - 1; j++) {
        f[j][0] = f[j][MESH_X - 2];
        f[j][MESH_X - 1] = f[j][1];
    }
    for (int i = 0; i < MESH_X - 1; i++) {
        f[0][i] = f[MESH_Y - 2][i];
        f[MESH_Y - 1][i] = f[1][i];
    }
}