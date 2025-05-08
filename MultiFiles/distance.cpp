#include "distance.h"
void distance(std::vector<Molecule>& mol, EM emf)
{
    double ro;
    double FLJ;
    double eps, sig;
    double kb = 1.38 * pow(10.0, -23);
    double kk = 8.98755 * pow(10.0, 9);
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            mol[i].AtomC[j].xa = mol[i].AtomC[j].x + mol[i].xc;
            mol[i].AtomC[j].ya = mol[i].AtomC[j].y + mol[i].yc;
            mol[i].AtomC[j].za = mol[i].AtomC[j].z + mol[i].zc;
            mol[i].AtomC[j].XX = 0;
            mol[i].AtomC[j].YY = 0;
            mol[i].AtomC[j].ZZ = 0;
        }
    }
    for (int i1 = 0; i1 < Molecule::numberofmolecules; i1++)
    {
       
        for (int j1 = 0; j1 < mol[i1].countofatoms; j1++)
        {
            for (int i2 = 0; i2 < Molecule::numberofmolecules; i2++)
            {
                if (i1 == i2) { continue; }
                for (int j2 = 0; j2 < mol[i2].countofatoms; j2++)
                {
                    ro = sqrt(pow(mol[i1].AtomC[j1].xa - mol[i2].AtomC[j2].xa, 2) + pow(mol[i1].AtomC[j1].ya - mol[i2].AtomC[j2].ya, 2) + pow(mol[i1].AtomC[j1].za - mol[i2].AtomC[j2].za, 2));
                    eps = sqrt(mol[i1].AtomC[j1].eps * mol[i2].AtomC[j2].eps) * kb;
                    sig = (mol[i1].AtomC[j1].sig + mol[i2].AtomC[j2].sig) / 2.0;
                    FLJ = 24.0 * eps / ro * pow(sig / ro, 6.0) * (2.0 * pow(sig / ro, 6) - 1.0); //–азмерность √игаЌьютон
                    FLJ += mol[i1].AtomC[j1].q * mol[i2].AtomC[j2].q / ro * kk;
                    mol[i1].AtomC[j1].XX += FLJ * (mol[i1].AtomC[j1].xa - mol[i2].AtomC[j2].xa) / ro; //размерность √игаЌьютон
                    mol[i1].AtomC[j1].YY += FLJ * (mol[i1].AtomC[j1].ya - mol[i2].AtomC[j2].ya) / ro; //размерность √игаЌьютон
                    mol[i1].AtomC[j1].ZZ += FLJ * (mol[i1].AtomC[j1].za - mol[i2].AtomC[j2].za) / ro; //размерность √игаЌьютон

                }
            }
        }
        
    }
    std::vector<double> Fxm, Fym, Fzm;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        
        Fxm.push_back(0.0); Fym.push_back(0.0); Fzm.push_back(0.0);
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if ((mol[i].AtomC[j].q == 0) || (emf.Bx == 0 && emf.By == 0 && emf.Bz == 0)) { continue; }
            Fxm[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].v * emf.Bz - mol[i].AtomC[j].w * emf.By);
            Fym[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].w * emf.Bx - mol[i].AtomC[j].u * emf.Bz);
            Fzm[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].u * emf.By - mol[i].AtomC[j].v * emf.Bx);
        }
        
        //cout<<Lxe[i]<<"\t"<<Lye[i]<<"\t"<<Lze[i]<<endl;
    }

    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        mol[i].ub = 0;
        mol[i].vb = 0;
        mol[i].wb = 0;
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            mol[i].ub += mol[i].AtomC[j].XX + mol[i].AtomC[j].q * emf.Ex * pow(10.0, -9.0)+ Fxm[i] * pow(10.0, -9.0);
            mol[i].vb += mol[i].AtomC[j].YY + mol[i].AtomC[j].q * emf.Ey * pow(10.0, -9.0)+ Fym[i] * pow(10.0, -9.0);
            mol[i].wb += mol[i].AtomC[j].ZZ + mol[i].AtomC[j].q * emf.Ez * pow(10.0, -9.0)+ Fzm[i] * pow(10.0, -9.0);
        }
        mol[i].ub /= mol[i].M;
        mol[i].vb /= mol[i].M;
        mol[i].wb /= mol[i].M;
        //cout<<mol[i].ub<<"\t"<<mol[i].vb<<"\t"<<mol[i].wb<<"\t"<<endl;
    }
    
}