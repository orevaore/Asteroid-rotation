#pragma once
#include<vector>
#include "Molecule.h"
#include "EM.h"
#include <fstream>
using namespace std;
void Epsil(std::vector<Molecule>& mol, double& Ep0, bool a, const EM& emf);