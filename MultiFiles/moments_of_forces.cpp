#include "moments_of_forces.h"
void moments_of_forces(std::vector<Molecule>& mol, EM const& emf)
{

    std::vector<double> Lx, Ly, Lz, Lxm, Lym, Lzm, Lxe, Lye, Lze,Lxm2,Lym2,Lzm2,Fxm,Fym,Fzm;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        Lx.push_back(0.0); Ly.push_back(0.0); Lz.push_back(0.0);
        Lxm.push_back(0.0); Lym.push_back(0.0); Lzm.push_back(0.0);
        Lxe.push_back(0.0); Lye.push_back(0.0); Lze.push_back(0.0);
        Lxm2.push_back(0.0); Lym2.push_back(0.0); Lzm2.push_back(0.0);
        Fxm.push_back(0.0); Fym.push_back(0.0); Fzm.push_back(0.0);
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            Lx[i] += mol[i].AtomC[j].y * mol[i].AtomC[j].ZZ - mol[i].AtomC[j].z * mol[i].AtomC[j].YY;
            Ly[i] += mol[i].AtomC[j].z * mol[i].AtomC[j].XX - mol[i].AtomC[j].x * mol[i].AtomC[j].ZZ;
            Lz[i] += mol[i].AtomC[j].x * mol[i].AtomC[j].YY - mol[i].AtomC[j].y * mol[i].AtomC[j].XX;
        }
    }
    //std::cout << Lx[2] << "\t" << Ly[2] << "\t" << Lz[2] << std::endl;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if ((mol[i].AtomC[j].q == 0.0) || (emf.Ex == 0.0 && emf.Ey == 0.0 && emf.Ez == 0.0)) { continue; }
            Lxe[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].y * emf.Ez - mol[i].AtomC[j].z * emf.Ey);
            Lye[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].z * emf.Ex - mol[i].AtomC[j].x * emf.Ez);
            Lze[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].x * emf.Ey - mol[i].AtomC[j].y * emf.Ex);
        }
        //cout<<Lxe[i]<<"\t"<<Lye[i]<<"\t"<<Lze[i]<<endl;
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if ((mol[i].AtomC[j].q == 0.0) || (emf.Bx == 0.0 && emf.By == 0.0 && emf.Bz == 0.0)) { continue; }
            Fxm[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].v * emf.Bz - mol[i].AtomC[j].w * emf.By);
            Fym[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].w * emf.Bx - mol[i].AtomC[j].u * emf.Bz);
            Fzm[i] += mol[i].AtomC[j].q * (mol[i].AtomC[j].u * emf.By - mol[i].AtomC[j].v * emf.Bx);
        }
        //cout<<Lxe[i]<<"\t"<<Lye[i]<<"\t"<<Lze[i]<<endl;
        
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if ((mol[i].AtomC[j].q == 0.0) || (emf.Bx == 0.0 && emf.By == 0.0 && emf.Bz == 0.0)) { continue; }
            Lxm2[i] +=  (mol[i].AtomC[j].y * Fzm[i] - mol[i].AtomC[j].z * Fym[i]);
            Lym2[i] +=  (mol[i].AtomC[j].z * Fxm[i] - mol[i].AtomC[j].x * Fzm[i]);
            Lzm2[i] +=  (mol[i].AtomC[j].x * Fym[i] - mol[i].AtomC[j].y * Fxm[i]);
        }
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if ((mol[i].AtomC[j].mux==0.0&& mol[i].AtomC[j].muy == 0.0&& mol[i].AtomC[j].muz == 0.0)) { continue; }
            Lxm[i] += mol[i].AtomC[j].muy * emf.Bz - mol[i].AtomC[j].muz * emf.By;
            Lym[i] += mol[i].AtomC[j].muz * emf.Bx - mol[i].AtomC[j].mux * emf.Bz;
            Lzm[i] += mol[i].AtomC[j].mux * emf.By - mol[i].AtomC[j].muy * emf.Bx;
            //std::cout << mol[i].AtomC[0].mux << "\t" << mol[i].AtomC[0].muy << "\t" << mol[i].AtomC[0].muz << std::endl;
        }
        //std::cout << Lxm[i] << "\t" << Lym[i] << "\t" << Lzm[i] << std::endl;
        
        //std::cout << emf.Bx << "\t" << emf.By << "\t" << emf.Bz << std::endl;
    }
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        mol[i].Lx = Lx[i] + Lxe[i] * pow(10.0, -9.0) + Lxm[i] +Lxm2[i] * pow(10.0, -9.0);
        mol[i].Ly = Ly[i] + Lye[i] * pow(10.0, -9.0) + Lym[i] +Lym2[i] * pow(10.0, -9.0);
        mol[i].Lz = Lz[i] + Lze[i] * pow(10.0, -9.0) + Lzm[i] +Lzm2[i] * pow(10.0, -9.0);
        // cout<<mol[i].Lx<<"\t"<<mol[i].Ly<<"\t"<<mol[i].Lz<<endl;
    }

}