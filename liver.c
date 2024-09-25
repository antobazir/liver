#include<stdio.h>
#include"structdef.h"
#include "initialize.h"
#include "diff_adv_reac.h"
#include "migrate.h"
/*code that leads other codes*/

Model Mod;
FILE *f_state;

void snap_print()
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
					Mod.M_Cell[k].G_cons,
					Mod.M_Cell[k].dG_cons,
					Mod.M_Cell[k].O_cons,
					Mod.M_Cell[k].dO_cons,
					Mod.M_Cell[k].G_diff,
					Mod.M_Cell[k].dG_diff,
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
						Mod.G.C[Mod.M_Cell[k].x_pos + Mod.M_Cell[k].y_pos*SZ1 + Mod.M_Cell[k].z_pos*SZ1*SZ2],
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


void main_loop(int Nmax)
{
    int i;

	snap_print();

	printf("main loop...\n");

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
            {
				migrate(i);
			}

            /*Kupffer cells*/
            if(Mod.M_Cell[i].Cell_type==4)
            {
				migrate(i);
			}

			/*divide(l+1);*/

			/*increases cell timer*/
			Mod.M_Cell[i].Timer = Mod.M_Cell[i].Timer+1;
        }

		

        snap_print();


        Mod.elapsed_mins++;
    } 
}


int main(int argc, char* argv[])
{
    float kG = 0.05;
    float kO = 0.05;

    f_state = fopen("liver_snap","w+");

    /*main code*/
    initialize(kG,kO);

    main_loop(10);

    fclose(f_state);


}