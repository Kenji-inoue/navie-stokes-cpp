#include "NavieStokes2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void initializeVelocity(Velocity2d& f, int meshX, int meshY, double lx, double ly) {
    const double PI=3.14159;
    const double K = 2 * PI;
    const double DX = lx / (meshX - 1);
    const double DY = ly / (meshY - 1);
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            const auto x = i * DX;
            const auto y = j * DY;
            f.u[j][i] = -1 * cos(K * x) * sin(K * y);
            f.v[j][i] = sin(K * x) * cos(K * y);
        }
    }
}

int main() {
    int meshX = 10;
    int meshY = meshX;
    double reynolds = 200;
    double lx = 1.0;
    double ly = lx;
    double deltaT = 0.02; 
    Velocity2d f;
    double omega = 1.0;
    double epsilon = 1e-7;
    double pRef = 1.0;
    MeshRange2d range = {1, meshX - 3, 1, meshY - 3};
    FieldUtil::setSize(f.u, meshX, meshY);
    FieldUtil::setSize(f.v, meshX, meshY);
    initializeVelocity(f, meshX, meshY, lx, ly);

    try {
        //NavieStokes2d navieStokes(meshX, meshY, reynolds, lx, ly, deltaT, f);
        const int interval = 1;
        const int maxIterations = 10;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(f.u, time, interval);
            FieldUtil::display(f.v, time, interval);
            //f = navieStokes.calculate();
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}