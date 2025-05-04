#include "NavieStokes2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void setInflowBoundaryCondition(Velocity2d& f, int meshX, int meshY) {
    for (int j = 0; j < meshY; j++) {
        f.u[j][0] = 0.98;
        f.u[j][1] = 0.98;
        f.v[j][0] = 0.02;
        f.v[j][1] = 0.02;
    }
    for (int i = 0; i < meshX; i++) {
        f.u[0][i] = 0.98;
        f.u[1][i] = 0.98;
        f.u[meshY - 1][i] = 0.98;
        f.u[meshY - 2][i] = 0.98;

        f.v[0][i] = 0.02;
        f.v[1][i] = 0.02;
        f.v[meshY - 1][i] = 0.02;
        f.v[meshY - 2][i] = 0.02;
    }
}

void defineObject(Object& object, int meshX, int meshY) {
    FieldUtil::InitializeFlagField(object.iu, meshX, meshY, ObjectFlag::fluid);
    FieldUtil::InitializeFlagField(object.ip, meshX, meshY, ObjectFlag::fluid);
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            if (28 <= j && j <= 36 && 28 <= i && i <= 36) {
                object.iu[j][i] = ObjectFlag::surface;
            }
            if (28 < j && j < 36 && 28 < i && i < 36) {
                object.iu[j][i] = ObjectFlag::inside;
            }
            if (28 <= j && j <= 35 && 28 <= i && i <= 35) {
                object.ip[j][i] = ObjectFlag::surface;
            }
            if (28 < j && j < 35 && 28 < i && i < 35) {
                object.ip[j][i] = ObjectFlag::inside;
            }
        }
    }
}

int main() {
    int meshX = 128 + 3;
    int meshY = 64 + 3;
    double reynolds = 200;
    double lx = 16.0;
    double dx = lx / (meshX - 3);
    double dy = dx;
    AnalysisResult result;
    double omega = 1.0;
    double epsilon = 1e-7;
    double pRef = 1.0;
    MeshRange2d range = {1, meshX - 3, 1, meshY - 3};
    FieldUtil::setSize(result.f.u, meshX, meshY);
    FieldUtil::setSize(result.f.v, meshX, meshY);
    FieldUtil::setSize(result.p, meshX, meshY);
    FieldUtil::setSize(result.s, meshX, meshY);
    FieldUtil::setSize(result.rot, meshX, meshY);
    setInflowBoundaryCondition(result.f, meshX, meshY);

    Object object;
    defineObject(object, meshX, meshY);

    try {
        NavieStokes2d solver(meshX, meshY, reynolds, dx, dy, 
                                  omega, epsilon, pRef, range, result, object);
        const int interval = 50;
        const int maxIterations = 500;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(result.rot, time, interval);
            result = solver.calculate();
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}