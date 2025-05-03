#include "NavieStokes2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void setInflowBoundaryCondition(Field2d& f, int meshX, int meshY) {
    for (int j = 0; j < meshY; j++) {
        f[j][0] = 1.0;
        f[j][1] = 1.0;
    }
}

int main() {
    int meshX = 16 + 3;
    int meshY = 8 + 3;
    double reynolds = 200;
    double lx = 1.0;
    double ly = lx;
    double dx = lx / (meshX - 3);
    double dy = ly / (meshY - 3);
    AnalysisResult result;
    double omega = 1.0;
    double epsilon = 1e-7;
    double pRef = 1.0;
    MeshRange2d range = {1, meshX - 3, 1, meshY - 3};
    FieldUtil::setSize(result.f.u, meshX, meshY);
    FieldUtil::setSize(result.f.v, meshX, meshY);
    FieldUtil::setSize(result.p, meshX, meshY);
    FieldUtil::setSize(result.s, meshX, meshY);
    setInflowBoundaryCondition(result.f.u, meshX, meshY);

    try {
        NavieStokes2d navieStokes(meshX, meshY, reynolds, dx, dy, 
                                  omega, epsilon, pRef, range, result);
        const int interval = 1;
        const int maxIterations = 2;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(result.f.u, time, interval);
            FieldUtil::display(result.f.v, time, interval);
            FieldUtil::display(result.p, time, interval);
            FieldUtil::display(result.s, time, interval);
            result = navieStokes.calculate();
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}