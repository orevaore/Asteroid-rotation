#include "Atom.h"
Atom::Atom()
{
    x = y = z = 0; /////////Временные координаты относительно центра масс
    xr1 = xr2 = xr3 = xr4 = yr1 = yr2 = yr3 = yr4 = zr1 = zr2 = zr3 = zr4 = 0; //координаты относительно центра масс
    xa = ya = za = 0; //координаты в абсолютном базисе
    ur1 = ur2 = ur3 = ur4 = vr1 = vr2 = vr3 = vr4 = wr1 = wr2 = wr3 = wr4 = ur = vr = wr=u=v=w = 0; //проекции скрорости атома
    XX = YY = ZZ = 0; //проекции сил
    mux = muy = muz = mux1 = muy1 = muz1 = muxb = muyb = muzb = 0;
    muxb1 = muxb2 = muxb3 = muxb4 = muyb1 = muyb2 = muyb3 = muyb4 = muzb1 = muzb2 = muzb3 = muzb4 = 0;
    m = 0.0035;
    eps = 5.0;
    sig = 0.34;
    q = 0;
    xa0 = ya0 = za0 = 0;
}