#pragma once
using namespace std;
#include <string>
#include "Atom.h"
enum xyz
{
    x, y, z
};

void ReadFile(const string directory, Atom* atom, xyz x);