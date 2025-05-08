#pragma once
#include "Atom.h"
#include <cmath>
class Atom
{
public:
    double x, y, z, m, eps, sig, q;
    double xr1, xr2, xr3, xr4, yr1, yr2, yr3, yr4, zr1, zr2, zr3, zr4; //координаты относительно центра масс
    double xa, ya, za; //координаты в абсолютном базисе
    double ur1, ur2, ur3, ur4, vr1, vr2, vr3, vr4, wr1, wr2, wr3, wr4, ur, vr, wr,u,v,w; //проекции скрорости атома
    double XX, YY, ZZ;//проекции сил
    double mux, muy, muz, mux1, muy1, muz1, muxb, muyb, muzb, muxb1, muxb2, muxb3, muxb4, muyb1, muyb2, muyb3, muyb4, muzb1, muzb2, muzb3, muzb4; //собственные магнитные моменты
    double xa0, ya0, za0;
    Atom();


};

