#include <string>
#include <stdexcept>
#include "NavieStokes2d.h"
#include "FieldUtil.h"

NavieStokes2d::NavieStokes2d(int meshX, int meshY, double reynolds, 
                        double dx, double dy, double omega, double epsilon, double pRef,
                        const MeshRange2d& range, AnalysisResult& result)
    : MESH_X(meshX), MESH_Y(meshY), REYNOLDS(reynolds),
      DX(dx), DY(dy), DELTA_T(0.2 * DX / 1.0),
      EPSILON(epsilon), P_REF(pRef), OMEGA(omega), MESH_RANGE(range), m_result(result),
      burgers_(meshX, meshY, reynolds, dx, dy, DELTA_T, result.f),
      poisson_(meshX, meshY, dx, dy, omega, epsilon, pRef, range)
{
    m_dp.resize(MESH_Y, std::vector<Value>(MESH_X, 0.0));
}

AnalysisResult NavieStokes2d::calculate() {
    m_result.f = calculateProvisionalVelocity(m_result.f, m_result.p);
    calculateDivergenceOfVelocity(m_result.s, m_result.f);
    const auto interval = poisson_.calculate(m_dp, m_result.s, 99999);
    modifyPressure(m_result.p, m_dp);
    modifyVelocity(m_result.f, m_dp);

    return m_result;
}

Velocity2d NavieStokes2d::calculateProvisionalVelocity(const Velocity2d& f, const Field2d& p) {
    Velocity2d f_next;
    FieldUtil::setSize(f_next.u, MESH_X, MESH_Y);
    FieldUtil::setSize(f_next.v, MESH_X, MESH_Y);

    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX+1; i <= MESH_RANGE.maxX; i++) {
            f_next.u[j][i] = f.u[j][i] + burgers_.calculateTerm(f.u, f, i, j) + calculatePressureTermX(p, i, j);
            f_next.v[j][i] = f.v[j][i] + burgers_.calculateTerm(f.v, f, i, j) + calculatePressureTermY(p, i, j);
        }
    }
    updateRunoffBoundaryCondition(f_next);
    return f_next;
}

Value NavieStokes2d::calculatePressureTermX(const Field2d& p, int i, int j) const {
    const Value frontPressure = (p[j][i] + p[j - 1][i]) / 2;
    const Value backPressure = (p[j][i - 1] + p[j - 1][i - 1]) / 2;
    return (frontPressure - backPressure) / DX;
}

Value NavieStokes2d::calculatePressureTermY(const Field2d& p, int i, int j) const {
    const Value frontPressure = (p[j][i] + p[j][i - 1]) / 2;
    const Value backPressure = (p[j - 1][i] + p[j - 1][i - 1]) / 2;
    return (frontPressure - backPressure) / DY;
}

void NavieStokes2d::calculateDivergenceOfVelocity(Field2d& s, const Velocity2d& f) {
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
            s[j][i] = ((f.u[j + 1][i + 1] - f.u[j + 1][i]) / DX + (f.u[j][i + 1] - f.u[j][i]) / DX) / 2 +
                      ((f.v[j + 1][i + 1] - f.v[j][i + 1]) / DY + (f.v[j + 1][i] - f.v[j][i]) / DY) / 2;
        }
    }
}

void NavieStokes2d::updateRunoffBoundaryCondition(Velocity2d& f) {
    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        f.u[j][MESH_X - 2] = f.u[j][MESH_X - 3];
        f.u[j][MESH_X - 1] = f.u[j][MESH_X - 3];
        f.v[j][MESH_X - 2] = f.v[j][MESH_X - 3];
        f.v[j][MESH_X - 1] = f.v[j][MESH_X - 3];
    }
}

void NavieStokes2d::modifyPressure(Field2d& p, Field2d& dp) {
    for (int j = 0; j < MESH_Y; j++) {
        for (int i = 0; i < MESH_X; i++) {
            p[j][i] += dp[j][i];
        }
    }

    // Runoff Boundary Condition
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        dp[j][MESH_X - 2] = p[j][MESH_X - 3] - p[j][MESH_X - 2];
        p[j][MESH_X - 2] = p[j][MESH_X - 3];
    }
    // Top and bottom edge
    for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
        p[0][i] = p[1][i];
        p[MESH_Y - 2][i] = p[MESH_Y - 3][i];
    }
}   

void NavieStokes2d::modifyVelocity(Velocity2d& f, const Field2d& dp) {
    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX+1; i <= MESH_RANGE.maxX; i++) {
            f.u[j][i] = f.u[j][i] - DELTA_T / 2 * ((dp[j][i] - dp[j][i - 1]) / DX + (dp[j - 1][i] - dp[j - 1][i - 1]) / DX);
            f.v[j][i] = f.v[j][i] - DELTA_T / 2 * ((dp[j][i] - dp[j - 1][i]) / DY + (dp[j][i - 1] - dp[j - 1][i - 1]) / DY);
        }
    }
}
