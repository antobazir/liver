#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structdef.h"

extern Model Mod;

void behav()
{
    int k;
    for(k=0;k<Mod.M_Tissue.N_Cell;k++)
    {
        /*hepatocyte only regulate blood sugar concentration*/
        if(Mod.M_Cell[k].Cell_type==1)
        {
            /*if hypoglaecimia*/
            if(Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]<)
            {

            }
        }

    }
}
