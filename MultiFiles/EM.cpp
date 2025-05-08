#include "EM.h"
EM::EM()
{
    Ex = 0;
    Ey = 0;
    Ez = 0;
    Bx = 0;
    By = 0;
    Bz = 0;
    E0 = 0; 
}
        void EM::parameters(double t)
        {
            double fr = 3.14 / 0.225      /2;         //         f/fr=1/4
            Ex = 0.0;// E0* cos(100 * t);
            Ey = 0;
            Ez = 0.0;
            Ez = 0.0;
            Bx = 0.0;
            //By = 20.0*cos(100*t);
           // Bz = 1.0*sin(fr*t)*1;
            Bz=0.0;
            //if (t > 1 / fr*4)
            //{
            //    Bz = 0.0;
            //}
            //if (t > 2.375)
            //{
            //    Bz = 0.0;
            //}
 /*           if (Bz<0)
            {
                Bz = 0.0;
            }*/
           // std::cout << Bz << std::endl;
        }