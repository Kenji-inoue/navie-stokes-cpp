#pragma once
#include "typedef.h"
#include "Advection2d.h"
#include "Diffusion2d.h"

class Burgers2d
{
public:
    Burgers2d(int meshX, int meshY, double reynolds, 
              double dx, double dy, double deltaT, Velocity2d f);
    ~Burgers2d() = default;

    Velocity2d calculate();
    Value calculateTerm(const Field2d& f, const Velocity2d& velocity, int i, int j) const;
    void validateTime();
    void updateBoundaryCondition(Field2d& f);
private:
    void updateVelocity(Velocity2d& f);
    const int MESH_X;
    const int MESH_Y;
    const double REYNOLDS;
    const double DX;
    const double DY;
    const double DELTA_T;
    Velocity2d m_f;
    Diffusion2d diffusion_;
    Advection2d advection_;
};
