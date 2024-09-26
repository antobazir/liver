#include<stdio.h>
#include"structdef.h"
#include "initialize.h"
#include "diff_adv_reac.h"
#include "behav.h"

/*code that leads other codes*/

Model Mod;
FILE *f_state,*f_G;

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

void print_slice_grid()
{

	int i,j;
	for(i=0;i<SZ1;i++)
	{
		for(j=0;j<SZ3;j++)
		{
			fprintf(f_G,"%3.3f ",Mod.G.C[i + SZ2/2*SZ1 + j*SZ1*SZ2]);
		}
		fprintf(f_G,"\n");
	}
	fprintf(f_G,"\n");
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

		
		/*the behavior function which changes the way cell behave*/
		behav();

		/*saves all the cells state variables in ascii*/
        snap_print();

		/*save central cross section along x-z plane*/
		print_slice_grid();

        Mod.elapsed_mins++;
    } 
}


int main(int argc, char* argv[])
{
    float kG = 0.01;
    float kO = 0.01;

    f_state = fopen("liver_snap","w+");
	f_G = fopen("G_snap","w+");

    /*main code*/
    initialize(kG,kO);

    main_loop(10);

    fclose(f_state);


}