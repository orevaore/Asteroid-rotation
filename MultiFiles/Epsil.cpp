#include "Epsil.h"
void Epsil(std::vector<Molecule>& mol, double& Ep0, bool a, const EM& emf)
{

    double ELJ = 0, alpha = 0, beta = 0, gamma = 0, MI = 0, Tb = 0, Ek = 0, omega = 0, Vel = 0;
    double epsil = 0;
    double EE = 0;
    double absVel = 0.0;
    ofstream F;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        omega = sqrt(pow(mol[i].omegax, 2) + pow(mol[i].omegay, 2) + pow(mol[i].omegaz, 2));

        if (omega == 0)
        {
            alpha = 0.0;
            beta = 0.0;
            gamma = 0.0;
        }
        else
        {
            alpha = mol[i].omegax / omega;
            beta = mol[i].omegay / omega;
            gamma = mol[i].omegaz / omega;
        }
        MI = mol[i].AA * alpha * alpha + mol[i].BB * beta * beta + mol[i].CC * gamma * gamma + 2 * mol[i].DD * beta * gamma + 2 * mol[i].EE * alpha * gamma + 2 * mol[i].FF * alpha * beta;
        Tb += 0.5 * MI * pow(omega, 2);
        Vel = sqrt(mol[i].u * mol[i].u + mol[i].v * mol[i].v + mol[i].w * mol[i].w);
        Ek += 0.5 * mol[i].M * pow(Vel, 2);
        if (i == Molecule::numberofmolecules) { continue; }
        double Veltemp = mol[i].u * mol[i].u + mol[i].v * mol[i].v + mol[i].w * mol[i].w;
        absVel += Veltemp;
        }
    absVel /= Molecule::numberofmolecules;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            mol[i].AtomC[j].xa = mol[i].AtomC[j].x + mol[i].xc;
            mol[i].AtomC[j].ya = mol[i].AtomC[j].y + mol[i].yc;
            mol[i].AtomC[j].za = mol[i].AtomC[j].z + mol[i].zc;
        }
    }
    if (a)
    {
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < mol[i].countofatoms; j++)
            {
                mol[i].AtomC[j].xa0 = mol[i].AtomC[j].x + mol[i].xc;
                mol[i].AtomC[j].ya0 = mol[i].AtomC[j].y + mol[i].yc;
                mol[i].AtomC[j].za0 = mol[i].AtomC[j].z + mol[i].zc;
            }
        }
    }


    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < mol[i].countofatoms; j++)
        {
            EE += mol[i].AtomC[j].q * (emf.Ex * (mol[i].AtomC[j].xa - mol[i].AtomC[j].xa0) + emf.Ey * (mol[i].AtomC[j].ya - mol[i].AtomC[j].ya0) + emf.Ez * (mol[i].AtomC[j].za - mol[i].AtomC[j].za0));
        }
    }
    double eps, sig, ro;
    for (int i1 = 0; i1 < Molecule::numberofmolecules; i1++)
    {
        for (int j1 = 0; j1 < mol[i1].countofatoms; j1++)
        {
            for (int i2 = 1; (i2 < Molecule::numberofmolecules); i2++)
            {
                if (i1 >= i2) { continue; }
                for (int j2 = 0; j2 < mol[i2].countofatoms; j2++)
                {
                    ro = sqrt(pow(mol[i1].AtomC[j1].xa - mol[i2].AtomC[j2].xa, 2) + pow(mol[i1].AtomC[j1].ya - mol[i2].AtomC[j2].ya, 2) + pow(mol[i1].AtomC[j1].za - mol[i2].AtomC[j2].za, 2));
                    eps = sqrt(mol[i1].AtomC[j1].eps * mol[i2].AtomC[j2].eps) * (1.38 * pow(10.0, -23));
                    sig = (mol[i1].AtomC[j1].sig + mol[i2].AtomC[j2].sig) / 2.0;
                    ELJ += 4.0 * eps * (pow((sig / ro), 12) - pow((sig / ro), 6));
                }
            }
        }

        if (a)
        {
            F.open("epsil.txt", ofstream::out | ofstream::trunc);
            F.close();
            F.open("alpha.txt", ofstream::out | ofstream::trunc);
            F.close();
            F.open("beta.txt", ofstream::out | ofstream::trunc);
            F.close();
            F.open("gamma.txt", ofstream::out | ofstream::trunc);
            F.close();
            F.open("absVel.txt", ofstream::out | ofstream::trunc);
            F.close();
            for (int i = 0; i < Molecule::numberofmolecules; i++)
            {

                Ep0 = Tb + Ek + ELJ;
            }
        }
    }
    epsil = (Tb + Ek + ELJ - Ep0) / (Ep0);
    //cout<<Ek<<"\t"<<Tb<<"\t"<<ELJ<<"\t"<<Ep0 << "\t" << EE <<endl;
    //cout<<" epslil =" <<epsil<<endl;

    if (a)
    {
        //epsil = 0;
    }
    F.open("epsil.txt", ios::app);
    F << epsil << "\n";
    F.close();
    F.open("alpha.txt", ios::app);
    F << alpha << "\n";
    F.close();
    F.open("beta.txt", ios::app);
    F << beta << "\n";
    F.close();
    F.open("gamma.txt", ios::app);
    F << gamma << "\n";
    F.close();
    F.open("absVel.txt", ios::app);
    F << absVel << "\n";
    F.close();
 
    F << epsil << "\n";
    F.close();
}