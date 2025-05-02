#include "Poisson2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>
#include <math.h>

void initializeSourceTerm(Field2d& s, int meshX, int meshY, double lx, double ly) {
    const double PI=3.14159;
    const double K_X = 2 * PI;
    const double K_Y = 2 * PI;
    const double DX = lx / (meshX - 1);
    const double DY = ly / (meshY - 1);
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            s[j][i] = -1 * (K_X * K_X + K_Y * K_Y) * sin(K_X * i * DX) * sin(K_Y * j * DY);
        }
    }
}

int main() {
    int meshX = 10;
    int meshY = meshX;
    double lx = 1.0;
    double ly = lx;
    double omega = 1.4;
    double epsilon = 1e-7;
    double pRef = 1.0;
    MeshRange2d range = {1, meshX - 2, 1, meshY - 2};
    Field2d s, p;
    FieldUtil::setSize(s, meshX, meshY);
    FieldUtil::setSize(p, meshX, meshY);
    initializeSourceTerm(s, meshX, meshY, lx, ly);

    try {
        Poisson2d poisson(meshX, meshY, lx, ly, omega, epsilon, pRef, range);
        const int interval = 10;
        const int maxIterations = 100;
        int iteration = 0;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(p, time, interval);
            iteration = poisson.calculate(p, s, maxIterations);
        }
        printf("Iteration: %d\n", iteration);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}