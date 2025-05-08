#include "urvrwr.h"
void urvrwr(std::vector<Molecule>& mol)
{
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            mol[i].AtomC[j].ur = mol[i].omegay * mol[i].AtomC[j].z - mol[i].omegaz * mol[i].AtomC[j].y;
            mol[i].AtomC[j].vr = mol[i].omegaz * mol[i].AtomC[j].x - mol[i].omegax * mol[i].AtomC[j].z;
            mol[i].AtomC[j].wr = mol[i].omegax * mol[i].AtomC[j].y - mol[i].omegay * mol[i].AtomC[j].x;
            mol[i].AtomC[j].u = mol[i].AtomC[j].ur + mol[i].u;
            mol[i].AtomC[j].v = mol[i].AtomC[j].vr + mol[i].v;
            mol[i].AtomC[j].w = mol[i].AtomC[j].wr + mol[i].w;
        }
    }
}