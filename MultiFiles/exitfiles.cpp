#include "exitfiles.h"
void exitfiles(vector<Molecule>& mol, bool a, double time0)
{
  //cout<<"here"<<endl;
    string txt = ".txt";
    string nit;
    string t = "t";
    ofstream f;
    vector<string> Vs;
    vector<ofstream> Vf;
    Vs.push_back("xcm");
    Vs.push_back("ycm");
    Vs.push_back("zcm");
    Vs.push_back("u");
    Vs.push_back("v");
    Vs.push_back("w");
    Vs.push_back("omegax");
    Vs.push_back("omegay");
    Vs.push_back("omegaz");
    Vs.push_back("x60");
    Vs.push_back("y60");
    Vs.push_back("z60");
    Vs.push_back("Velocity");
    Vs.push_back("omegamod");
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < Vs.size(); j++)
        {
            Vf.push_back(ofstream());
        }
    }
    f.open(t + txt, ios::app);
    f << time0 << "\n";
    
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        stringstream ss;
        ss << i;
        nit = ss.str();
        Vf[0].open(Vs[0] + nit + txt, ios::app);
        Vf[0] << mol[i].xc << "\n";
        Vf[0].close();
        Vf[1].open(Vs[1] + nit + txt, ios::app);
        Vf[1] << mol[i].yc << "\n";
        Vf[1].close();
        Vf[2].open(Vs[2] + nit + txt, ios::app);
        Vf[2] << mol[i].zc << "\n";
        Vf[2].close();
        Vf[3].open(Vs[3] + nit + txt, ios::app);
        Vf[3] << mol[i].u << "\n";
        Vf[3].close();
        Vf[4].open(Vs[4] + nit + txt, ios::app);
        Vf[4] << mol[i].v << "\n";
        Vf[4].close();
        Vf[5].open(Vs[5] + nit + txt, ios::app);
        Vf[5] << mol[i].w << "\n";
        Vf[5].close();
        Vf[6].open(Vs[6] + nit + txt, ios::app);
        Vf[6] << mol[i].omegax << "\n";
        Vf[6].close();
        Vf[7].open(Vs[7] + nit + txt, ios::app);
        Vf[7] << mol[i].omegay << "\n";
        Vf[7].close();
        Vf[8].open(Vs[8] + nit + txt, ios::app);
        Vf[8] << mol[i].omegaz << "\n";
        Vf[8].close();
        Vf[9].open(Vs[9] + nit + txt, ios::app);
        Vf[9] << mol[i].AtomC[17].x << "\n";
        Vf[9].close();
        Vf[10].open(Vs[10] + nit + txt, ios::app);
        Vf[10] << mol[i].AtomC[17].y << "\n";
        Vf[10].close();
        Vf[11].open(Vs[11] + nit + txt, ios::app);
        Vf[11] << mol[i].AtomC[17].z << "\n";
        Vf[11].close();
        Vf[11].open(Vs[12] + nit + txt, ios::app);
        Vf[11] << sqrt(mol[i].AtomC[0].u* mol[i].AtomC[0].u + mol[i].AtomC[0].v* mol[i].AtomC[0].v + mol[i].AtomC[0].w* mol[i].AtomC[0].w) << "\n";
        Vf[11].close();
        Vf[11].open(Vs[13] + nit + txt, ios::app);
        Vf[11] << sqrt(mol[i].omegax * mol[i].omegax + mol[i].omegay * mol[i].omegay + mol[i].omegaz * mol[i].omegaz) << "\n";
        Vf[11].close();

    }

}