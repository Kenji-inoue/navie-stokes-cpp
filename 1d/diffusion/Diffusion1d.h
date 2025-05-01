#pragma once

class Diffusion1d
{
public:
    Diffusion1d();
    ~Diffusion1d() = default;
    void validateTime();
    void initializeField();

    void displayProcess(int time);
    void calculate();
    void updateField();
    void simulate();

};