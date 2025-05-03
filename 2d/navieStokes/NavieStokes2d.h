#pragma once
#include "typedef.h"
#include "Burgers2d.h"
#include "Poisson2d.h"

class NavieStokes2d
{
public:
    NavieStokes2d(int meshX, int meshY, double reynolds, 
                double lx, double ly, double deltaT, Velocity2d f,
                double omega, double epsilon, double pRef,
                const MeshRange2d& range, Field2d& p, Field2d& s);
    ~NavieStokes2d() = default;
    void caclulate();

    
private:
    Velocity2d calculateProvisionalVelocity(const Velocity2d& f, const Field2d& p);
    Value calculatePressureTermX(const Field2d& p, int i, int j) const;
    Value calculatePressureTermY(const Field2d& p, int i, int j) const;
    void calculateDivergenceOfVelocity(Field2d& s, const Velocity2d& f);
    void updateRunoffBoundaryCondition(Velocity2d& f);
    void modifyPressure(Field2d& p, Field2d& dp);
    void modifyVelocity(Velocity2d& f, const Field2d& dp);
    const int MESH_X;
    const int MESH_Y;
    const double REYNOLDS;
    const double DX;
    const double DY;
    const double DELTA_T;
    const double EPSILON;
    const double P_REF;
    const double OMEGA;
    const MeshRange2d MESH_RANGE;
    Velocity2d m_f;
    Field2d m_p;
    Field2d m_dp;
    Field2d m_s;
    Burgers2d burgers_;
    Poisson2d poisson_;
};
