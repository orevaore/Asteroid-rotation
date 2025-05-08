#include "tenzors.h"

void tenzors(std::vector<Molecule>& mol, bool a)
{
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {

        mol[i].AA = 0;
        mol[i].BB = 0;
        mol[i].CC = 0;
        mol[i].DD = 0;
        mol[i].EE = 0;
        mol[i].FF = 0;
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            mol[i].AA += (pow(mol[i].AtomC[j].y, 2) + pow(mol[i].AtomC[j].z, 2)) * mol[i].AtomC[j].m;
            mol[i].BB += (pow(mol[i].AtomC[j].z, 2) + pow(mol[i].AtomC[j].x, 2)) * mol[i].AtomC[j].m;
            mol[i].CC += (pow(mol[i].AtomC[j].x, 2) + pow(mol[i].AtomC[j].y, 2)) * mol[i].AtomC[j].m;
            //         mol[i].DD-=mol[i].AtomC[j].y*mol[i].AtomC[j].z*mol[i].AtomC[j].m;
            //         mol[i].EE-=mol[i].AtomC[j].z*mol[i].AtomC[j].x*mol[i].AtomC[j].m;
            //         mol[i].FF-=mol[i].AtomC[j].x*mol[i].AtomC[j].y*mol[i].AtomC[j].m;
            mol[i].DD += (pow(mol[i].AtomC[j].y - mol[i].AtomC[j].z, 2) - (pow(mol[i].AtomC[j].y, 2) + pow(mol[i].AtomC[j].z, 2))) * mol[i].AtomC[j].m / 2;
            mol[i].EE += (pow(mol[i].AtomC[j].z - mol[i].AtomC[j].x, 2) - (pow(mol[i].AtomC[j].z, 2) + pow(mol[i].AtomC[j].x, 2))) * mol[i].AtomC[j].m / 2;
            mol[i].FF += (pow(mol[i].AtomC[j].x - mol[i].AtomC[j].y, 2) - (pow(mol[i].AtomC[j].x, 2) + pow(mol[i].AtomC[j].y, 2))) * mol[i].AtomC[j].m / 2;
            //  cout<<mol[i].AtomC[j].x-mol[i].AtomC[j].xr2<<"  "<<mol[i].AtomC[j].y-mol[i].AtomC[j].yr2<<"  "<<mol[i].AtomC[j].z-mol[i].AtomC[j].zr2<<endl;
            // cout<<mol[i].AtomC[j].x-mol[i].AtomC[j].xr1<<"  "<<mol[i].AtomC[j].y-mol[i].AtomC[j].yr1<<"  "<<mol[i].AtomC[j].z-mol[i].AtomC[j].zr1<<endl;
         // cout<<mol[i].AtomC[j].x-mol[i].AtomC[j].xr3<<"  "<<mol[i].AtomC[j].y-mol[i].AtomC[j].yr3<<"  "<<mol[i].AtomC[j].z-mol[i].AtomC[j].zr3<<endl;
        }
     //std::cout<<"A="<<mol[i].AA<<"  B="<<mol[i].BB<<"  C="<<mol[i].CC<<std::endl;
        //    cout<<"D="<<mol[i].DD<<"  E="<<mol[i].EE<<"  F="<<mol[i].FF<<endl;
    }
    std::vector<double> omznam, omchx, omchy, omchz;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        omznam.push_back(0); omchx.push_back(0); omchy.push_back(0); omchz.push_back(0);
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        if (a) { continue; }
        omznam[i] = mol[i].AA * mol[i].DD * mol[i].DD - 2 * mol[i].DD * mol[i].EE * mol[i].FF + mol[i].BB * mol[i].EE * mol[i].EE + mol[i].CC * mol[i].FF * mol[i].FF - mol[i].AA * mol[i].BB * mol[i].CC;
        omchx[i] = mol[i].DD * mol[i].DD * mol[i].Kx - mol[i].BB * mol[i].CC * mol[i].Kx + mol[i].BB * mol[i].EE * mol[i].Kz + mol[i].CC * mol[i].FF * mol[i].Ky - mol[i].DD * mol[i].EE * mol[i].Ky - mol[i].DD * mol[i].FF * mol[i].Kz;
        omchy[i] = mol[i].EE * mol[i].EE * mol[i].Ky - mol[i].AA * mol[i].CC * mol[i].Ky + mol[i].AA * mol[i].DD * mol[i].Kz + mol[i].CC * mol[i].FF * mol[i].Kx - mol[i].DD * mol[i].EE * mol[i].Kx - mol[i].EE * mol[i].FF * mol[i].Kz;
        omchz[i] = mol[i].FF * mol[i].FF * mol[i].Kz - mol[i].AA * mol[i].BB * mol[i].Kz + mol[i].AA * mol[i].DD * mol[i].Ky + mol[i].BB * mol[i].EE * mol[i].Kx - mol[i].DD * mol[i].FF * mol[i].Kx - mol[i].EE * mol[i].FF * mol[i].Ky;

        mol[i].omegax = omchx[i] / omznam[i];
        mol[i].omegay = omchy[i] / omznam[i];
        mol[i].omegaz = omchz[i] / omznam[i];
        if (omznam[i] == 0)
        {
            mol[i].omegax = mol[i].omegay = mol[i].omegaz = 0.0;
        }
        //     cout<<"tenzors"<<endl;
        //    std::cout<<mol[i].omegax<<"  "<<mol[i].omegay<<"  "<<mol[i].omegaz<<"  "<<std::endl;
        //     cout<<omchx[i]<<"  "<<omchy[i]<<"  "<<omchz[i]<<"  "<<endl;
        //     cout<<mol[i].AA<<"  "<<mol[i].BB<<"  "<<mol[i].CC<<endl;
    }
}