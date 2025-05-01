#pragma once
#include "typedef.h"

class Advection2d
{
public:
    Advection2d(int meshX, int meshY, double constU, double constV, double deltaX, double deltaY, double deltaT);
    ~Advection2d() = default;

    Field2d calculate(const Field2d& f);
private:
    void validateTime();
    void updateBoundaryCondition(Field2d& f);
    const int MESH_X;
    const int MESH_Y;
    const double CONST_U;
    const double CONST_V;
    const double DELTA_X;
    const double DELTA_Y;
    const double DELTA_T;
};
