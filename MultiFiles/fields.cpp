#include "fields.h"
void fields(std::vector<Molecule>& mol, double dt, bool a)
{
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            if (mol[i].AtomC[j].mux || mol[i].AtomC[j].muy || mol[i].AtomC[j].muz)
            {
                if (a)
                {
                    mol[i].AtomC[j].mux = mol[i].AtomC[j].mux1 + dt * mol[i].AtomC[j].muxb;
                    mol[i].AtomC[j].muy = mol[i].AtomC[j].muy1 + dt * mol[i].AtomC[j].muyb;
                    mol[i].AtomC[j].muz = mol[i].AtomC[j].muz1 + dt * mol[i].AtomC[j].muzb;
                }
                mol[i].AtomC[j].muxb = mol[i].omegay * mol[i].AtomC[j].muz - mol[i].omegaz * mol[i].AtomC[j].muy;
                mol[i].AtomC[j].muyb = mol[i].omegaz * mol[i].AtomC[j].mux - mol[i].omegax * mol[i].AtomC[j].muz;
                mol[i].AtomC[j].muzb = mol[i].omegax * mol[i].AtomC[j].muy - mol[i].omegay * mol[i].AtomC[j].mux;
                //cout<<mol[i].AtomC[j].mux<<"  "<<mol[i].AtomC[j].muy<<"  "<<mol[i].AtomC[j].muz<<endl;
                //cout<<sqrt(pow(mol[i].AtomC[j].mux,2)+pow(mol[i].AtomC[j].muy,2)+pow(mol[i].AtomC[j].muz,2))<<endl;
                //cout<<mol[i].AtomC[j].muxb<<"  "<<mol[i].AtomC[j].muyb<<"  "<<mol[i].AtomC[j].muzb<<endl;
            }

        }
    }
}