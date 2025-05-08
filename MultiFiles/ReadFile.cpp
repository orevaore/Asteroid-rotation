
#include "ReadFile.h"
#include <fstream>


void ReadFile(const string directory, Atom* atom, xyz x)
{
    int i = 0;
    ifstream File;
    File.open(directory);
    while (true)
        if (!File.eof())
        {
            if (x == xyz::x)
            {
                File >> atom[i].x;
            }
            else if (x == xyz::y)
            {
                File >> atom[i].y;
            }
            else if (x == xyz::z)
            {
                File >> atom[i].z;
            }
            i++;
        }
        else
        {
            break;
        }

    File.close();
}