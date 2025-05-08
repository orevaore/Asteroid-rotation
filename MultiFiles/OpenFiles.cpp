#include "OpenFiles.h"

int OpenFiles()
{
    
    string txt = ".txt";
    string nit;
    string t = "t";

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
    ofstream f;
    f.open(t+txt, ofstream::out | ofstream::trunc);
    f.close();
    int countoffiles = 0;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        stringstream ss;
        ss << i;
        nit = ss.str();
        for (int j = 0; j < Vs.size(); j++)
        {
            Vf[countoffiles].open(Vs[j] + nit + txt, ofstream::out | ofstream::trunc);
            Vf[countoffiles].close();
            countoffiles++;
        }
    }


   
    return (Vs.size());
}
