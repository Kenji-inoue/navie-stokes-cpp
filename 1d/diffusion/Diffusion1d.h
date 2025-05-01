#pragma once
#include <vector>

class Diffusion1d
{
public:
    Diffusion1d(int meshX, double constA, double deltaX, double deltaT);
    ~Diffusion1d() = default;
    void validateTime();
    void initializeField();

    void displayProcess(int time);
    void calculate();
    void updateField();
    void simulate();

private:
    const int MESH_X;
    const double CONST_A;
    const double DELTA_X;
    const double DELTA_T;

    std::vector<double> f;
    std::vector<double> f_next;

};
