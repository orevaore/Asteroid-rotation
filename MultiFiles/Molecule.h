#pragma once
#include "Atom.h"
#include <vector>
class Molecule
{
public:
    static int numberofmolecules;
    double M;
    double xc, yc, zc, u, v, w, u1, v1, w1, omegax, omegay, omegaz, omega, Vel;
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4; //Положения центра масс в 4 промежуточных позициях
    double AA, BB, CC, DD, EE, FF; //компоненты тензора инерции
    double Kx, Ky, Kz, Kx1, Ky1, Kz1; //проекции кинетического момента на оси абсолютного базиса(текущий и первый на каждом полном шаге
    double Lx1, Lx2, Lx3, Lx4, Ly1, Ly2, Ly3, Ly4, Lz1, Lz2, Lz3, Lz4, Lx, Ly, Lz;//Проекции моментов сил(полные)
    double ub1, ub2, ub3, ub4, vb1, vb2, vb3, vb4, wb1, wb2, wb3, wb4, ub, vb, wb; //ускрорения
    double xb1, xb2, xb3, xb4, yb1, yb2, yb3, yb4, zb1, zb2, zb3, zb4;
    double alpha, beta, gamma, Ek, Tb, MI, ELJ; //Величины, в основном для проверки баланса энергии
    enum xyz
    {
        x, y, z
    };
    Atom* AtomC;
    int countofatoms;


    Molecule(int countofatoms);
    void cmtozero();
    void rotate(double degree, xyz x);
    void changem(int indexofatom, double newm);


};

