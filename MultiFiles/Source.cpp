//Инструкция
//чтобы добавить молекулу существующего типа, необходимо изменить добавить numberofsamemol[i]=1; Задать начальные данные для молекулы;
//Для добавление нового типа молекул необходимо 1. увеличить numberdifferencemol=1; 2. Затем указать новое numberofsamemol[i]=1; 3. Указать число атомов в молекуле numberofatoms[i]=61; 
// 4. Указать новый путь directories[0][0]="xC.txt"; directories[0][1]="yC.txt"; directories[0][2]="zC.txt"; 5. Начальные данные.
//массу атомов необходимо менять через метод Molecule void changem(int indexofatom, double newm)
//добавить private поля для безопасности. Поместить туда массу...
#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <chrono>
#include "ReadFile.h"
#include "Atom.h"
#include "Em.h"
#include "Molecule.h"
#include "distance.h"
#include "tenzors.h"
#include "urvrwr.h"
#include "fields.h"
#include "moments_of_forces.h"
#include "Epsil.h"
#include "OpenFiles.h"
#include "exitfiles.h"
//void ft (double &u, double &v, double &w, double dt, int steps);
using namespace std;
const int stepsnumber = 1000000*1, numberdifferencemol = 1, eachstep = 2000;
const double dt = pow(10.0, -6.0);
int main()
{
    auto start = chrono::high_resolution_clock::now();

    EM emf;

    bool del = true;
    int numberofsamemol[numberdifferencemol] = {}; /////don't touch count type of molecules
    string directories[numberdifferencemol][3] = {};////don't touch
    int numberofatoms[numberdifferencemol] = {};   ////don't touch
    //////////initial+numberdifferencemol
    numberofsamemol[0] = 1;
    numberofatoms[0] = 100;
    for (int i = 0; i < numberdifferencemol; i++)
    {
        directories[i][0] = "xC" + to_string(i) + ".txt";
        directories[i][1] = "yC" + to_string(i) + ".txt";
        directories[i][2] = "zC" + to_string(i) + ".txt";
    }



    //////////initial+numberdifferencemol
    vector<Molecule> vMol;
    int tempiter = 0;
    for (int i = 0; i < numberdifferencemol; i++)
    {
        for (int j = 0; j < numberofsamemol[i]; j++)
        {
            vMol.push_back(Molecule(numberofatoms[i]));
            ReadFile(directories[i][0], vMol[tempiter].AtomC, x);
            ReadFile(directories[i][1], vMol[tempiter].AtomC, y);
            ReadFile(directories[i][2], vMol[tempiter].AtomC, z);
            tempiter++;
        }
    }

    vector<ofstream> Vf0;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        vMol[i].cmtozero();
    }
    int numberoffiles = OpenFiles();
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < numberoffiles; j++)
        {
            Vf0.push_back(ofstream());
        }
    }

    double E0;
    //initial data






    vMol[0].omegax = 0.0;
    vMol[0].omegay = 100.0;
    vMol[0].omegaz = 0.0;
    vMol[0].u = 0.0;
    vMol[0].v = 0.0;
    vMol[0].w = 0.0;
    vMol[0].rotate(75, Molecule::z);
 
    
      //vMol[0].changem(23, pow(10.0, -26) * 9.27);
     //vMol[0].changem(29, pow(10.0, -26) * 9.27);
    //vMol[0].changem(34, pow(10.0, -26) * 9.27);
   // vMol[0].AtomC[0].q = 1.6 * pow(10.0, -19.0);
  //vMol[0].AtomC[0].muz = 9.274 * pow(10.0, -24.0) * 3.23;



// 
   // vMol[1].AtomC[240].q = 1.6 * pow(10.0, -19.0);
    //vMol[1].AtomC[0].mux=9.274*pow(10.0,-24.0)*3.23; ///Магнетон Бора - элементарный магнитный момент. Магнитный момент железа ~3.23 Магнетона Бора
 /*   vMol[0].xc = 1.1;
  
     vMol[2].omegax=500.0;*/
     //  vMol[1].omegay=10.0;
     //  vMol[1].omegaz=1.0;
    

     // vMol[1].xc=2.0;



      //initial data

    tenzors(vMol, true);
    //cout<<endl<<vMol[0].AA<<"\t"<<vMol[0].BB<<"\t"<<vMol[0].CC<<endl<<vMol[0].DD<<"\t"<<vMol[0].EE<<"\t"<<vMol[0].FF<<"\t"<<endl;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        vMol[i].Kx = vMol[i].AA * vMol[i].omegax + vMol[i].FF * vMol[i].omegay + vMol[i].EE * vMol[i].omegaz;
        vMol[i].Ky = vMol[i].FF * vMol[i].omegax + vMol[i].BB * vMol[i].omegay + vMol[i].DD * vMol[i].omegaz;
        vMol[i].Kz = vMol[i].EE * vMol[i].omegax + vMol[i].DD * vMol[i].omegay + vMol[i].CC * vMol[i].omegaz;
        //     cout<<"0\t"<<vMol[i].Kx<<" "<<vMol[i].Ky<<" "<<vMol[i].Kz<<endl;
    }
    /*   for (int i = 0; i < Molecule::numberofmolecules; i++)
       {
           for (int j = 0; j < vMol[i].countofatoms; j++)
           {
               vMol[i].AtomC[j].xa0 = vMol[i].AtomC[j].x + vMol[i].xc;
               vMol[i].AtomC[j].ya0 = vMol[i].AtomC[j].y + vMol[i].yc;
               vMol[i].AtomC[j].za0 = vMol[i].AtomC[j].z + vMol[i].zc;
           }
       }*/
    setlocale(LC_ALL, "rus");
    cout << "Всего типов молекул = " << numberdifferencemol << endl;
    cout << "Количество одинаковых молекул \t";
    for (int i = 0; i < numberdifferencemol; i++)
    {
        cout << numberofsamemol[i] << "\t";
    }
    cout << endl << "Начальные скорости\t";
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        cout << vMol[i].u << "\t" << vMol[i].v << "\t" << vMol[i].w << "\t" << endl;
    }
    cout << endl << "Начальные вращения\t";
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        cout << vMol[i].omegax << "\t" << vMol[i].omegay << "\t" << vMol[i].omegaz << "\t" << endl;
    }
    cout << endl << "Величины полей\t";




    cout << endl << "Наличие зарядов:" << endl;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < vMol[i].countofatoms; j++)
        {
            if (vMol[i].AtomC[j].q)
            {
                cout << "Величина заряда=" << vMol[i].AtomC[j].q << endl;
                cout << "Заряд находится в " << i << " молекуле. Индекс атома " << j << endl;
            }
        }
    }
    cout << endl << "Собственные магнитные моменты:" << endl;
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        for (int j = 0; j < vMol[i].countofatoms; j++)
        {
            if (vMol[i].AtomC[j].mux || vMol[i].AtomC[j].muy || vMol[i].AtomC[j].muz)
            {
                cout << "Величины собственных магнитных моментов=" << sqrt(pow(vMol[i].AtomC[j].mux, 2) + pow(vMol[i].AtomC[j].muy, 2) + pow(vMol[i].AtomC[j].muz, 2)) << endl;
                cout << "собственным магнитным моментом обладает атом в " << i << " молекуле. Индекс атома " << j << endl;
            }
        }
    }
    cout << "Проверьте пожалуйста emf.parameters" << endl;
    cout << "Число шагов по времени = " << stepsnumber << endl;
    cout << "Величина шага = " << dt << endl;
    cout << endl;

    double time = 0;
    for (int steps = 0; steps < stepsnumber; steps++) ///////////main cycle 
    {
        if (!(steps % eachstep))
        {
            cout << steps * 100 / stepsnumber << "%" << "\r";
        }

        time = double(steps) * dt;                  ////////////calculate time every step

        emf.parameters(time);
        urvrwr(vMol);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].Kx1 = vMol[i].Kx;
            vMol[i].Ky1 = vMol[i].Ky;
            vMol[i].Kz1 = vMol[i].Kz;
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].ur1 = vMol[i].AtomC[j].ur;
                vMol[i].AtomC[j].vr1 = vMol[i].AtomC[j].vr;
                vMol[i].AtomC[j].wr1 = vMol[i].AtomC[j].wr;
                //  cout<<vMol[i].AtomC[j].ur1<<" "<<vMol[i].AtomC[j].vr1<<" "<< vMol[i].AtomC[j].wr1<<endl;
            }
        }
        distance(vMol, emf);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].ub1 = vMol[i].ub;
            vMol[i].vb1 = vMol[i].vb;
            vMol[i].wb1 = vMol[i].wb;
            //    cout<<vMol[i].ub1<<" "<<vMol[i].vb1<<" "<<vMol[i].wb1<<endl;
        }
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].mux1 = vMol[i].AtomC[j].mux;
                vMol[i].AtomC[j].muy1 = vMol[i].AtomC[j].muy;
                vMol[i].AtomC[j].muz1 = vMol[i].AtomC[j].muz;
            }
        }

        fields(vMol, dt / 2, false);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].muxb1 = vMol[i].AtomC[j].muxb;
                vMol[i].AtomC[j].muyb1 = vMol[i].AtomC[j].muyb;
                vMol[i].AtomC[j].muzb1 = vMol[i].AtomC[j].muzb;
            }
        }
        moments_of_forces(vMol, emf);

        //RungeKutta2
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].u1 = vMol[i].u;//?
            vMol[i].v1 = vMol[i].v;//?
            vMol[i].w1 = vMol[i].w;//?
            vMol[i].xb1 = vMol[i].u;//?
            vMol[i].yb1 = vMol[i].v;//?
            vMol[i].zb1 = vMol[i].w;//?
            vMol[i].x1 = vMol[i].xc;//?
            vMol[i].y1 = vMol[i].yc;//?
            vMol[i].z1 = vMol[i].zc;//?=

            vMol[i].x2 = vMol[i].x1 + dt / 2.0 * vMol[i].u1;
            vMol[i].y2 = vMol[i].y1 + dt / 2.0 * vMol[i].v1;
            vMol[i].z2 = vMol[i].z1 + dt / 2.0 * vMol[i].w1;

            vMol[i].u = vMol[i].u1 + dt / 2.0 * vMol[i].ub1;
            vMol[i].v = vMol[i].v1 + dt / 2.0 * vMol[i].vb1;
            vMol[i].w = vMol[i].w1 + dt / 2.0 * vMol[i].wb1;
            vMol[i].xb2 = vMol[i].u;
            vMol[i].yb2 = vMol[i].v;
            vMol[i].zb2 = vMol[i].w;
            vMol[i].Lx1 = vMol[i].Lx;
            vMol[i].Ly1 = vMol[i].Ly;
            vMol[i].Lz1 = vMol[i].Lz;
            vMol[i].Kx = vMol[i].Lx1 * dt / 2.0 + vMol[i].Kx1;
            vMol[i].Ky = vMol[i].Ly1 * dt / 2.0 + vMol[i].Ky1;
            vMol[i].Kz = vMol[i].Lz1 * dt / 2.0 + vMol[i].Kz1;
            //         cout<<"1\t"<<vMol[i].Kx<<"\t"<<vMol[i].Ky<<"\t"<<vMol[i].Kz<<endl;
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].xr1 = vMol[i].AtomC[j].x;
                vMol[i].AtomC[j].yr1 = vMol[i].AtomC[j].y;
                vMol[i].AtomC[j].zr1 = vMol[i].AtomC[j].z;
                vMol[i].AtomC[j].xr2 = vMol[i].AtomC[j].xr1 + dt / 2.0 * vMol[i].AtomC[j].ur;
                vMol[i].AtomC[j].yr2 = vMol[i].AtomC[j].yr1 + dt / 2.0 * vMol[i].AtomC[j].vr;
                vMol[i].AtomC[j].zr2 = vMol[i].AtomC[j].zr1 + dt / 2.0 * vMol[i].AtomC[j].wr;
                vMol[i].AtomC[j].x = vMol[i].AtomC[j].xr2;
                vMol[i].AtomC[j].y = vMol[i].AtomC[j].yr2;
                vMol[i].AtomC[j].z = vMol[i].AtomC[j].zr2;
                // cout<<vMol[i].AtomC[j].xr2<<"\t"<< vMol[i].AtomC[j].yr2<<"\t"<< vMol[i].AtomC[j].zr2<<endl;
            }
            vMol[i].xc = vMol[i].x2;
            vMol[i].yc = vMol[i].y2;
            vMol[i].zc = vMol[i].z2;
        }
        tenzors(vMol, false);
        urvrwr(vMol);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].ur2 = vMol[i].AtomC[j].ur;
                vMol[i].AtomC[j].vr2 = vMol[i].AtomC[j].vr;
                vMol[i].AtomC[j].wr2 = vMol[i].AtomC[j].wr;
                //cout<<vMol[i].AtomC[j].ur2<<"\t"<< vMol[i].AtomC[j].vr2<<"\t"<< vMol[i].AtomC[j].wr2<<endl;
            }
        }
        emf.parameters(time + dt / 2);
        distance(vMol, emf);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].ub2 = vMol[i].ub;
            vMol[i].vb2 = vMol[i].vb;
            vMol[i].wb2 = vMol[i].wb;
            //  cout<<vMol[i].ub2<<" "<<vMol[i].vb2<<" "<<vMol[i].wb2<<endl;
        }
        fields(vMol, dt / 2, true);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].muxb2 = vMol[i].AtomC[j].muxb;
                vMol[i].AtomC[j].muyb2 = vMol[i].AtomC[j].muyb;
                vMol[i].AtomC[j].muzb2 = vMol[i].AtomC[j].muzb;
            }
        }
        moments_of_forces(vMol, emf);
        //RungeKutta3
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].x3 = vMol[i].x1 + dt / 2.0 * vMol[i].xb2;
            vMol[i].y3 = vMol[i].y1 + dt / 2.0 * vMol[i].yb2;
            vMol[i].z3 = vMol[i].z1 + dt / 2.0 * vMol[i].zb2;
            vMol[i].u = vMol[i].u1 + dt / 2.0 * vMol[i].ub2;
            vMol[i].v = vMol[i].v1 + dt / 2.0 * vMol[i].vb2;
            vMol[i].w = vMol[i].w1 + dt / 2.0 * vMol[i].wb2;
            vMol[i].xb3 = vMol[i].u;
            vMol[i].yb3 = vMol[i].v;
            vMol[i].zb3 = vMol[i].w;
            vMol[i].Lx2 = vMol[i].Lx;
            vMol[i].Ly2 = vMol[i].Ly;
            vMol[i].Lz2 = vMol[i].Lz;
            vMol[i].Kx = vMol[i].Lx2 * dt / 2.0 + vMol[i].Kx1;
            vMol[i].Ky = vMol[i].Ly2 * dt / 2.0 + vMol[i].Ky1;
            vMol[i].Kz = vMol[i].Lz2 * dt / 2.0 + vMol[i].Kz1;
            //cout<<vMol[i].u<<" "<<vMol[i].v<<" "<<vMol[i].w<<endl;
   //         cout<<"2\t" <<vMol[i].Kx<<"\t"<< vMol[i].Ky<<"\t"<< vMol[i].Kz<<endl;
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].xr3 = vMol[i].AtomC[j].xr1 + dt / 2.0 * vMol[i].AtomC[j].ur2;
                vMol[i].AtomC[j].yr3 = vMol[i].AtomC[j].yr1 + dt / 2.0 * vMol[i].AtomC[j].vr2;
                vMol[i].AtomC[j].zr3 = vMol[i].AtomC[j].zr1 + dt / 2.0 * vMol[i].AtomC[j].wr2;
                vMol[i].AtomC[j].x = vMol[i].AtomC[j].xr3;
                vMol[i].AtomC[j].y = vMol[i].AtomC[j].yr3;
                vMol[i].AtomC[j].z = vMol[i].AtomC[j].zr3;
                //   cout<<vMol[i].AtomC[j].xr3<<"  "<<vMol[i].AtomC[j].yr3<<"  "<<vMol[i].AtomC[j].zr3<<endl;
            }
            vMol[i].xc = vMol[i].x3;
            vMol[i].yc = vMol[i].y3;
            vMol[i].zc = vMol[i].z3;
        }
        tenzors(vMol, false);
        urvrwr(vMol);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].ur3 = vMol[i].AtomC[j].ur;
                vMol[i].AtomC[j].vr3 = vMol[i].AtomC[j].vr;
                vMol[i].AtomC[j].wr3 = vMol[i].AtomC[j].wr;
                //cout<<vMol[i].AtomC[j].ur2<<"\t"<< vMol[i].AtomC[j].vr2<<"\t"<< vMol[i].AtomC[j].wr2<<endl;
            }
        }
        emf.parameters(time + dt / 2);
        distance(vMol, emf);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].ub3 = vMol[i].ub;
            vMol[i].vb3 = vMol[i].vb;
            vMol[i].wb3 = vMol[i].wb;
            //  cout<<vMol[i].ub3<<" "<<vMol[i].vb3<<" "<<vMol[i].wb3<<endl;
        }
        fields(vMol, dt / 2, true);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].muxb3 = vMol[i].AtomC[j].muxb;
                vMol[i].AtomC[j].muyb3 = vMol[i].AtomC[j].muyb;
                vMol[i].AtomC[j].muzb3 = vMol[i].AtomC[j].muzb;
            }
        }
        moments_of_forces(vMol, emf);
        //RungeKutta4
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].x4 = vMol[i].x1 + dt * vMol[i].xb3;
            vMol[i].y4 = vMol[i].y1 + dt * vMol[i].yb3;
            vMol[i].z4 = vMol[i].z1 + dt * vMol[i].zb3;
            vMol[i].u = vMol[i].u1 + dt * vMol[i].ub3;
            vMol[i].v = vMol[i].v1 + dt * vMol[i].vb3;
            vMol[i].w = vMol[i].w1 + dt * vMol[i].wb3;
            vMol[i].xb4 = vMol[i].u;
            vMol[i].yb4 = vMol[i].v;
            vMol[i].zb4 = vMol[i].w;
            vMol[i].Lx3 = vMol[i].Lx;
            vMol[i].Ly3 = vMol[i].Ly;
            vMol[i].Lz3 = vMol[i].Lz;
            vMol[i].Kx = vMol[i].Lx3 * dt + vMol[i].Kx1;
            vMol[i].Ky = vMol[i].Ly3 * dt + vMol[i].Ky1;
            vMol[i].Kz = vMol[i].Lz3 * dt + vMol[i].Kz1;
            //  cout<<vMol[i].u<<" "<<vMol[i].v<<" "<<vMol[i].w<<endl;
    //         cout<<"3\t"<< vMol[i].Kx<<"\t"<<vMol[i].Ky<<"\t"<< vMol[i].Kz<<endl;
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].xr4 = vMol[i].AtomC[j].xr1 + dt * vMol[i].AtomC[j].ur3;
                vMol[i].AtomC[j].yr4 = vMol[i].AtomC[j].yr1 + dt * vMol[i].AtomC[j].vr3;
                vMol[i].AtomC[j].zr4 = vMol[i].AtomC[j].zr1 + dt * vMol[i].AtomC[j].wr3;
                vMol[i].AtomC[j].x = vMol[i].AtomC[j].xr4;
                vMol[i].AtomC[j].y = vMol[i].AtomC[j].yr4;
                vMol[i].AtomC[j].z = vMol[i].AtomC[j].zr4;
            }
            //   cout<<vMol[i].omegax<<"  "<<vMol[i].omegay<<"  "<<vMol[i].omegaz<<endl;
            vMol[i].xc = vMol[i].x4;
            vMol[i].yc = vMol[i].y4;
            vMol[i].zc = vMol[i].z4;
            //cout<< vMol[i].Kx<<" "<< vMol[i].Ky<<" "<< vMol[i].Kz<<endl;
        }
        tenzors(vMol, false);
        urvrwr(vMol);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].ur4 = vMol[i].AtomC[j].ur;
                vMol[i].AtomC[j].vr4 = vMol[i].AtomC[j].vr;
                vMol[i].AtomC[j].wr4 = vMol[i].AtomC[j].wr;
                //cout<<vMol[i].AtomC[j].ur2<<"\t"<< vMol[i].AtomC[j].vr2<<"\t"<< vMol[i].AtomC[j].wr2<<endl;
            }
        }
        emf.parameters(time + dt);
        distance(vMol, emf);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].ub4 = vMol[i].ub;
            vMol[i].vb4 = vMol[i].vb;
            vMol[i].wb4 = vMol[i].wb;
            //   cout<<vMol[i].ub4<<" "<<vMol[i].vb4<<" "<<vMol[i].wb4<<endl;
        }
        fields(vMol, dt, true);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].muxb4 = vMol[i].AtomC[j].muxb;
                vMol[i].AtomC[j].muyb4 = vMol[i].AtomC[j].muyb;
                vMol[i].AtomC[j].muzb4 = vMol[i].AtomC[j].muzb;
            }
        }
        moments_of_forces(vMol, emf);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].Lx4 = vMol[i].Lx;
            vMol[i].Ly4 = vMol[i].Ly;
            vMol[i].Lz4 = vMol[i].Lz;
        }
        //RungeKutta full value
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            vMol[i].u = vMol[i].u1 + dt / 6.0 * (vMol[i].ub1 + 2.0 * vMol[i].ub2 + 2.0 * vMol[i].ub3 + vMol[i].ub4);
            vMol[i].v = vMol[i].v1 + dt / 6.0 * (vMol[i].vb1 + 2.0 * vMol[i].vb2 + 2.0 * vMol[i].vb3 + vMol[i].vb4);
            vMol[i].w = vMol[i].w1 + dt / 6.0 * (vMol[i].wb1 + 2.0 * vMol[i].wb2 + 2.0 * vMol[i].wb3 + vMol[i].wb4);
            vMol[i].xc = vMol[i].x1 + dt / 6.0 * (vMol[i].xb1 + 2.0 * vMol[i].xb2 + 2.0 * vMol[i].xb3 + vMol[i].xb4);
            vMol[i].yc = vMol[i].y1 + dt / 6.0 * (vMol[i].yb1 + 2.0 * vMol[i].yb2 + 2.0 * vMol[i].yb3 + vMol[i].yb4);
            vMol[i].zc = vMol[i].z1 + dt / 6.0 * (vMol[i].zb1 + 2.0 * vMol[i].zb2 + 2.0 * vMol[i].zb3 + vMol[i].zb4);
            vMol[i].Kx = vMol[i].Kx1 + dt / 6.0 * (vMol[i].Lx1 + 2.0 * vMol[i].Lx2 + 2.0 * vMol[i].Lx3 + vMol[i].Lx4);
            vMol[i].Ky = vMol[i].Ky1 + dt / 6.0 * (vMol[i].Ly1 + 2.0 * vMol[i].Ly2 + 2.0 * vMol[i].Ly3 + vMol[i].Ly4);
            vMol[i].Kz = vMol[i].Kz1 + dt / 6.0 * (vMol[i].Lz1 + 2.0 * vMol[i].Lz2 + 2.0 * vMol[i].Lz3 + vMol[i].Lz4);
            //ЗАКОН ОТРАЖЕНИЯ
            double polkuba = 5.0;
            //if (vMol[i].xc> polkuba && vMol[i].u>0 || vMol[i].xc < -polkuba && vMol[i].u < 0)
            //{ 
            //    vMol[i].u = -vMol[i].u;
            //}
            //if (vMol[i].yc > polkuba && vMol[i].v > 0 || vMol[i].yc < -polkuba && vMol[i].v < 0)
            //{
            //    vMol[i].v = -vMol[i].v;
            //}
            //if (vMol[i].zc > polkuba && vMol[i].w > 0 || vMol[i].zc < -polkuba && vMol[i].w < 0)
            //{
            //    vMol[i].w = -vMol[i].w;
            //}
               //ЗАКОН ОТРАЖЕНИЯ
            for (int j = 0; j < vMol[i].countofatoms; j++)
            {
                vMol[i].AtomC[j].x = vMol[i].AtomC[j].xr1 + dt / 6 * (vMol[i].AtomC[j].ur1 + 2 * vMol[i].AtomC[j].ur2 + 2 * vMol[i].AtomC[j].ur3 + vMol[i].AtomC[j].ur4);
                vMol[i].AtomC[j].y = vMol[i].AtomC[j].yr1 + dt / 6 * (vMol[i].AtomC[j].vr1 + 2 * vMol[i].AtomC[j].vr2 + 2 * vMol[i].AtomC[j].vr3 + vMol[i].AtomC[j].vr4);
                vMol[i].AtomC[j].z = vMol[i].AtomC[j].zr1 + dt / 6 * (vMol[i].AtomC[j].wr1 + 2 * vMol[i].AtomC[j].wr2 + 2 * vMol[i].AtomC[j].wr3 + vMol[i].AtomC[j].wr4);
                vMol[i].AtomC[j].mux = vMol[i].AtomC[j].mux1 + dt / 6 * (vMol[i].AtomC[j].muxb1 + 2 * vMol[i].AtomC[j].muxb2 + 2 * vMol[i].AtomC[j].muxb3 + vMol[i].AtomC[j].muxb4);
                vMol[i].AtomC[j].muy = vMol[i].AtomC[j].muy1 + dt / 6 * (vMol[i].AtomC[j].muyb1 + 2 * vMol[i].AtomC[j].muyb2 + 2 * vMol[i].AtomC[j].muyb3 + vMol[i].AtomC[j].muyb4);
                vMol[i].AtomC[j].muz = vMol[i].AtomC[j].muz1 + dt / 6 * (vMol[i].AtomC[j].muzb1 + 2 * vMol[i].AtomC[j].muzb2 + 2 * vMol[i].AtomC[j].muzb3 + vMol[i].AtomC[j].muzb4);


            }
            //  cout<<vMol[i].AtomC[60].mux<<"  "<<vMol[i].AtomC[60].muy<< "  "<<vMol[i].AtomC[60].muz<<endl;
              /*
          cout<<vMol[i].AtomC[60].muxb1<<"  "<<vMol[i].AtomC[60].muyb1<< "  "<<vMol[i].AtomC[60].muzb1<<endl;
          cout<<vMol[i].AtomC[60].muxb2<<"  "<<vMol[i].AtomC[60].muyb2<< "  "<<vMol[i].AtomC[60].muzb2<<endl;
          cout<<vMol[i].AtomC[60].muxb3<<"  "<<vMol[i].AtomC[60].muyb3<< "  "<<vMol[i].AtomC[60].muzb3<<endl;
          cout<<vMol[i].AtomC[60].muxb4<<"  "<<vMol[i].AtomC[60].muyb4<< "  "<<vMol[i].AtomC[60].muzb4<<endl;*/
        }
        tenzors(vMol, false);
        for (int i = 0; i < Molecule::numberofmolecules; i++)
        {
            //cout<<vMol[i].omegax<<" "<<vMol[i].omegay<<" "<<vMol[i].omegaz<<endl;
        }
        if (!(steps % eachstep))
        {
            exitfiles(vMol, del, time);
            Epsil(vMol, E0, del, emf);
        }

        del = false;
    }
    ofstream F;
    F.open("VelAll.txt", ofstream::out | ofstream::trunc);
    F.close();
    F.open("VelAll.txt", ios::app);
    for (int i = 0; i < Molecule::numberofmolecules; i++)
    {
        F << sqrt(vMol[i].u * vMol[i].u + vMol[i].v * vMol[i].v + vMol[i].w * vMol[i].w) << "\n";
    }
    F.close();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<float>duration = end - start;
    cout << "duration = " << duration.count() << endl;
    return 0;
}

