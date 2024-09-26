#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structdef.h"
#include "migrate.h"

extern Model Mod;

void behav()
{
    int k;
    for(k=0;k<Mod.M_Tissue.N_Cell;k++)
    {
        /*for heaptocytes*/
        if(Mod.M_Cell[k].Cell_type==1)
        {
            /*if hypoglaecimia && cell not already responding*/
            if(
                (Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]<3.7)
                &&(Mod.M_Cell[k].Chg_timer==0)
                ) /*3.7 mM threshold source : PMID 22783644 */
            {
                Mod.M_Cell[k].Chg_timer= 5*Mod.M_Tissue.reac_time; /*15 mn delay + 60mn to completely switch to sugar delivery*/
                Mod.M_Cell[k].dG_cons= ((-Mod.kG) - Mod.M_Cell[k].G_cons)/4*Mod.M_Tissue.reac_time;  
            }

            /*if homeoglaecimia && cell not already responding*/
            if(
                (Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]<4.2)
                &&(Mod.M_Cell[k].Chg_timer==0)
                ) /*3.7 mM threshold source : PMID 22783644 */
            {
                Mod.M_Cell[k].Chg_timer= 5*Mod.M_Tissue.reac_time; /*15 mn delay + 60mn to completely switch to sugar delivery*/  
                Mod.M_Cell[k].dG_cons= (Mod.kG - Mod.M_Cell[k].G_cons)/4*Mod.M_Tissue.reac_time;  
            }
        }

        /*for all cells currently responding */
        if(Mod.M_Cell[k].Chg_timer>0)
        {
            /*decrement timer*/
            Mod.M_Cell[k].Chg_timer = Mod.M_Cell[k].Chg_timer-1;

            /*if the 15 minutes delay for response is finished start modfiying glucose consumpyion*/
            if((Mod.M_Cell[k].Cell_type==1)&&
                (Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]<3.7)
                &&(Mod.M_Cell[k].Chg_timer<=60)) 
            {
                /*change the glucose consumption gradually and update the grid*/
                Mod.M_Cell[k].G_cons = Mod.M_Cell[k].G_cons + Mod.M_Cell[k].dG_cons;
                Mod.G.kC[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2] = Mod.M_Cell[k].G_cons;
            }

                        /*if the 15 minutes delay for response is finished start modfiying glucose consumpyion*/
            if((Mod.M_Cell[k].Cell_type==1)&&
                (Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]>=4.5)
                &&(Mod.M_Cell[k].Chg_timer<=60)) 
            {
                /*change the glucose consumption gradually and update the grid*/
                Mod.M_Cell[k].G_cons = Mod.M_Cell[k].G_cons + Mod.M_Cell[k].dG_cons;
                Mod.G.kC[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2] = Mod.M_Cell[k].G_cons;
            }

        }

        /*move stellate cells*/
        if(Mod.M_Cell[k].Cell_type==3)
        {
            migrate(k);
        }

        /*Kupffer cells*/
        if(Mod.M_Cell[k].Cell_type==4)
        {
            migrate(k);
        } 

        /*increases cell timer*/
        Mod.M_Cell[k].Timer = Mod.M_Cell[k].Timer+1;   

    }
}
