#include "Molecule.h"
#include <iostream>
int Molecule::numberofmolecules = 0;
Molecule::Molecule(int countofatoms = 0)
{
    this->countofatoms = countofatoms;
    numberofmolecules++;
    xc = yc = zc = 0; //временные координаты центра масс
    x1 = x2 = x3 = x4 = y1 = y2 = y3 = y4 = z1 = z2 = z3 = z4 = 0;//Положения центра масс в 4 промежуточных позициях
    M = 0;
    u = v = w = Vel = u1 = v1 = w1 = 0;
    ub1 = ub2 = ub3 = ub4 = vb1 = vb2 = vb3 = vb4 = wb1 = wb2 = wb3 = wb4 = ub = vb = wb = 0;
    omegax = omegay = omegaz = omega = 0;
    AA = BB = CC = DD = EE = FF = Kx = Ky = Kz = Kx1 = Ky1 = Kz1 = 0;
    Lx1 = Lx2 = Lx3 = Lx4 = Ly1 = Ly2 = Ly3 = Ly4 = Lz1 = Lz2 = Lz3 = Lz4 = Lx = Ly = Lz = 0;
    alpha = beta = gamma = Ek = Tb = MI = ELJ = 0;
    xb1 = xb2 = xb3 = xb4 = yb1 = yb2 = yb3 = yb4 = zb1 = zb2 = zb3 = zb4 = 0;
    AtomC = new Atom[countofatoms]; 
    for (int i = 0; i < countofatoms; i++)
    {
        M += AtomC[i].m;
    }
}
void Molecule::cmtozero()
{
    double deltacmx = 0, deltacmy = 0, deltacmz = 0;
    for (int i = 0; i < this->countofatoms; i++)
    {
        deltacmx += this->AtomC[i].x * this->AtomC[i].m;
        deltacmy += this->AtomC[i].y * this->AtomC[i].m;
        deltacmz += this->AtomC[i].z * this->AtomC[i].m;
    }
    
    deltacmx /= this->M;
    deltacmy /= this->M;
    deltacmz /= this->M;
    std::cout << "dcmx = " << deltacmx << "\tdcmy = " << deltacmy << "\tdcmz = " << deltacmz << std::endl;

        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].x -= deltacmx;
        }

        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].y -= deltacmy;
        }
        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].z -= deltacmz;
        }

}

void Molecule::rotate(double degree, xyz x)
{
    double radian = degree / 180 * (atan(1.0) * 4.0);
    std::vector<double> xtemp, ytemp, ztemp;
    for (int i = 0; i < countofatoms; i++)
    {
        xtemp.push_back(0); ytemp.push_back(0); ztemp.push_back(0);
    }
    if (x == xyz::x)
    {
        for (int i = 0; i < this->countofatoms; i++)
        {
            ytemp[i] = this->AtomC[i].y * cos(radian) - this->AtomC[i].z * sin(radian);
            ztemp[i] = this->AtomC[i].y * sin(radian) + this->AtomC[i].z * cos(radian);
        }
        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].y = ytemp[i];
            this->AtomC[i].z = ztemp[i];
        }

    }
    else if (x == xyz::y)
    {
        for (int i = 0; i < this->countofatoms; i++)
        {
            xtemp[i] = this->AtomC[i].x * cos(radian) + this->AtomC[i].z * sin(radian);
            ztemp[i] = -this->AtomC[i].x * sin(radian) + this->AtomC[i].z * cos(radian);
        }
        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].x = xtemp[i];
            this->AtomC[i].z = ztemp[i];
        }
    }
    else if (x == xyz::z)
    {
        for (int i = 0; i < this->countofatoms; i++)
        {
            xtemp[i] = this->AtomC[i].x * cos(radian) - this->AtomC[i].y * sin(radian);
            ytemp[i] = this->AtomC[i].x * sin(radian) + this->AtomC[i].y * cos(radian);
        }
        for (int i = 0; i < this->countofatoms; i++)
        {
            this->AtomC[i].x = xtemp[i];
            this->AtomC[i].y = ytemp[i];
        }
    }
    cmtozero();
}

void Molecule::changem(int indexofatom, double newm)
{
    this->M -= this->AtomC[indexofatom].m - newm;
    this->AtomC[indexofatom].m = newm;
    cmtozero();
}
