#include <string>
#include <stdexcept>
#include "NavieStokes2d.h"
#include "FieldUtil.h"

NavieStokes2d::NavieStokes2d(int meshX, int meshY, double reynolds, 
                     double lx, double ly, double deltaT, Velocity2d f) 
    : MESH_X(meshX), MESH_Y(meshY), REYNOLDS(reynolds), 
      DX(lx / (meshX - 1)), DY(ly / (meshY - 1)), DELTA_T(deltaT), m_f(f),
      diffusion_(meshX, meshY, 1/reynolds, DX, DY, deltaT), 
      advection_(meshX, meshY, DX, DY, deltaT)
{
    validateTime();
}

void NavieStokes2d::validateTime() {
    diffusion_.validateTime();

    const auto maxU = FieldUtil::findMax(m_f.u);
    const auto maxV = FieldUtil::findMax(m_f.v);
    advection_.validateTime(maxU, maxV);
}

Velocity2d NavieStokes2d::calculate() {
    Velocity2d f_next;
    FieldUtil::setSize(f_next.u, MESH_X, MESH_Y);
    FieldUtil::setSize(f_next.v, MESH_X, MESH_Y);

    for (int j = 1; j <= MESH_Y - 2; j++) {
        for (int i = 1; i <= MESH_X - 2; i++) {
            f_next.u[j][i] = m_f.u[j][i] + calculateTerm(m_f.u, m_f, i, j);
            f_next.v[j][i] = m_f.v[j][i] + calculateTerm(m_f.v, m_f, i, j);
        }
    }
    updateBoundaryCondition(f_next.u);
    updateBoundaryCondition(f_next.v);
    updateVelocity(f_next);
    return m_f;
}

Value NavieStokes2d::calculateTerm(const Field2d& f, const Velocity2d& velocity, int i, int j) const {
    return advection_.calculateVelocity(f, velocity.u, velocity.v, i, j) + diffusion_.calculateTerm(f, i, j);
}

void NavieStokes2d::updateBoundaryCondition(Field2d& f) {
    advection_.updateBoundaryCondition(f);
}

void NavieStokes2d::updateVelocity(Velocity2d& f) {
    std::swap(m_f, f);
}
