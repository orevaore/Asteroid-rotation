#pragma once
#include <cmath>
#include <iostream>
class EM
{
public:
    double Ex, Ey, Ez, Bx, By, Bz, E0;
    EM();
    void  parameters(double t);
};

