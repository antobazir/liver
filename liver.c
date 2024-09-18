#include<stdio.h>
#include"structdef.h"
#include "initialize.h"
#include "diff_adv_reac.h"
/*code that leads other codes*/

Model Mod;
FILE *f_state;

void main_loop(int Nmax)
{
    int i;

    while(Mod.elapsed_mins<Nmax)
    {
        /*running diffusion advection for one minute*/
        diff_adv_reac();

        for(i=0;i<Mod.M_Tissue.N_Cell;i++)
        {
            /*hepatocyte*/
            if(Mod.M_Cell[i].Cell_type==1)
            {}

            /*lsec*/
            if(Mod.M_Cell[i].Cell_type==2)
            {}
            
            /*stellate cells*/
            if(Mod.M_Cell[i].Cell_type==3)
            {}

            /*Kupffer cells*/
            if(Mod.M_Cell[i].Cell_type==4)
            {}
        }

        log_print();

        Mod.elapsed_mins++;
    } 
}

void log_print()
{
    int k;

    for(k=0;k<Mod.M_Tissue.N_Cell;k++)
				{
					/*printf("%d\n",k);*/
					fprintf(f_state,"%d %d %d %d %f %f %f %f %f %f %f %f %d %d %d %d %d",
					Mod.M_Cell[k].Cell_index,
                    Mod.M_Cell[k].Timer,
					Mod.M_Cell[k].Cycle_dur,
					Mod.M_Cell[k].Chg_timer,
					Mod.M_Cell[k].S_cons,
					Mod.M_Cell[k].dS_cons,
					Mod.M_Cell[k].O_cons,
					Mod.M_Cell[k].dO_cons,
					Mod.M_Cell[k].S_diff,
					Mod.M_Cell[k].dS_diff,
					Mod.M_Cell[k].O_diff,
					Mod.M_Cell[k].dO_diff,
					Mod.M_Cell[k].x_pos,
					Mod.M_Cell[k].y_pos,
					Mod.M_Cell[k].z_pos,
					Mod.M_Cell[k].parent_idx,
					Mod.M_Cell[k].Cell_type);

					if((Mod.M_Cell[k].x_pos!=-1)&&(Mod.M_Cell[k].y_pos!=-1)&&(Mod.M_Cell[k].z_pos!=-1))
					{
						fprintf(f_state," %1.4f %1.4f %d %d\n",
						Mod.S.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2],
						Mod.O.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2],
						Mod.M_Tissue.Grid[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2],
						Mod.M_Tissue.LD[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2]
						);
						
					}
					else
					{
						fprintf(f_state," %1.4f %1.4f %d %d\n",
						0.0,
						0.0,
						0,
						0
						);
						
					}
				
				}
				fprintf(f_state,"\n");
}

int main(int argc, char* argv[])
{
    float kS = 1.0;
    float kO = 1.0;

    f_state = fopen("liver_snap","w+");

    /*main code*/
    initialize(kS,kO);

    main_loop(1);

    fclose(f_state);


}