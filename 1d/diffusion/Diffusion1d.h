#pragma once
#include "typedef.h"

class Diffusion1d
{
public:
    Diffusion1d(int meshX, double constA, double deltaX, double deltaT);
    ~Diffusion1d() = default;

    void displayProcess(int time, int interval, const Field1d& f);
    Field1d calculate(const Field1d& f);

private:
    void initializeField();
    void validateTime();
    const int MESH_X;
    const double CONST_A;
    const double DELTA_X;
    const double DELTA_T;

    Field1d f_next;
};
