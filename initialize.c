#include<stdio.h>
#include<stdlib.h>
#include"structdef.h"/*the header containing model types*/


extern Model Mod;/*line that link this invokes the model instance called in liver.c*/

void place_hepatocyte(int row, int col, int stck)
{
    /*increasing cell number*/
    Mod.M_Tissue.N_Cell++;

    /*placing things on Grids*/
    Mod.M_Tissue.Grid[row + col*SZ1 +stck*SZ1*SZ2] =  Mod.M_Tissue.N_Cell;/*first line*/
    Mod.G.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.G.DC_tiss;
    Mod.G.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kG;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    /*updating the cell array*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 1; /*hepatocyte cell type=1*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_cons = Mod.kG;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_diff = Mod.G.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = Mod.kO;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_diff = Mod.O.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].parent_idx = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Chg_timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cycle_dur = 1500;
}



void place_Disse(int row, int col, int stck)
{
    Mod.G.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.G.DC_mat;
    Mod.G.kC[row + col*SZ1 +stck*SZ1*SZ2] = 0.0;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_mat;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = 0.0;
}

void place_lsec(int row, int col, int stck)
{
    /*increasing cell number*/
    Mod.M_Tissue.N_Cell++;

    /*placing things on Grids*/
    Mod.M_Tissue.Grid[row + col*SZ1 +stck*SZ1*SZ2] =  Mod.M_Tissue.N_Cell;
    Mod.G.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.G.DC_endo;
    Mod.G.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kG;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_endo;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 2; /*lse cell type=2*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_cons = Mod.kG;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_diff = Mod.G.DC_endo;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = Mod.kO;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_diff = Mod.O.DC_endo;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].parent_idx = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Chg_timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cycle_dur = 1500;
}

void place_stellate(int row, int col, int stck)
{
    /*increasing cell number*/
    Mod.M_Tissue.N_Cell++;

    /*placing things on Grids*/
    Mod.M_Tissue.Grid[row + col*SZ1 +stck*SZ1*SZ2] =  Mod.M_Tissue.N_Cell;
    Mod.G.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.G.DC_tiss;
    Mod.G.kC[row + col*SZ1 +stck*SZ1*SZ2] = 2*Mod.kG;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = 2*Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 3; /*stellate cell type=3*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_cons = 2*Mod.kG;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_diff = Mod.G.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = 2*Mod.kO;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_diff = Mod.O.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].parent_idx = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Chg_timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cycle_dur = 1500;
}

void place_Kupffer(int row, int col, int stck)
{
    /*increasing cell number*/
    Mod.M_Tissue.N_Cell++;

    /*placing things on Grids*/
    Mod.M_Tissue.Grid[row + col*SZ1 +stck*SZ1*SZ2] =  Mod.M_Tissue.N_Cell;
    Mod.G.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.G.DC_tiss;
    Mod.G.kC[row + col*SZ1 +stck*SZ1*SZ2] = 2*Mod.kG;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 4; /*Kupffer cell type=4*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_cons = 2*Mod.kG;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].G_diff = Mod.G.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dG_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = 2*Mod.kO;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_diff = Mod.O.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].parent_idx = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Chg_timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cycle_dur = 1500;
}

void initialize(float kG, float kO)
{
    printf("initializing... \n");

    Mod.kO = kO;
	Mod.kG = kG;
    FILE* file_field;

    int i=0,j=0,k=0,l=0;

    /*setting constants: units | concentration : mM | length : µm | */
    Mod.M_Tissue.N_Cell =0;
    Mod.G.DC_med = 700.0;
    Mod.G.DC_mat = 500.0;
    Mod.G.DC_tiss = 120.0;
    Mod.G.DC_endo = 50.0;
    Mod.O.DC_med = 3500.0;
    Mod.O.DC_mat = 3000.0;
    Mod.O.DC_tiss = 2000.0;
    Mod.O.DC_endo = 2000.0;
    Mod.G.C_ext = 3.5;
    Mod.O.C_ext = 0.07;
    
    Mod.M_Tissue.reac_time = 15*60; /*source 10.1677/joe.0.1260109*/


    /*allocating memory for the nutrient grids*/
    Mod.G.C = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.G.DCm = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.G.kC = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.G.fl_x = malloc(SZ2*SZ3*sizeof(double));

    /*initializing variables*/
    Mod.G.dx = 2.0;
    Mod.G.dt =0.25*Mod.G.dx*Mod.G.dx/Mod.G.DC_med;

    Mod.O.C = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.DCm = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.kC = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.fl_x = malloc(SZ2*SZ3*sizeof(double));
    
    Mod.O.dx = 2.0;
    Mod.O.dt =0.25*Mod.O.dx*Mod.O.dx/Mod.O.DC_med;
    
    /*Grid to follow cells*/
    Mod.M_Tissue.Grid = malloc(SZ1*SZ2*SZ3*sizeof(int));
	Mod.M_Tissue.LD = malloc(SZ1*SZ2*SZ3*sizeof(int));
	Mod.M_Tissue.state_mat = malloc(SZ1*SZ2*SZ3*sizeof(unsigned char));
	Mod.M_Cell = malloc(300000*sizeof(Cell));


    /*read the velocity field in the file generated by the GNU/Octave script : square_poiseuille_flow.m*/
    file_field = fopen("velocity_field","r");
    fread(Mod.G.fl_x,sizeof(double),SZ2*SZ3,file_field);

    for(k=0;k<SZ3;k++)
    {
        for(j=0;j<SZ2;j++)
        {
              for(i=0;i<SZ1;i++)
            {
                Mod.G.DCm[i + j*SZ1 + k*SZ1*SZ2] = Mod.G.DC_med;
                Mod.O.DCm[i + j*SZ1 + k*SZ1*SZ2] = Mod.O.DC_med;
                Mod.G.C[i + j*SZ1 + k*SZ1*SZ2] = Mod.G.C_ext;
                Mod.O.C[i + j*SZ1 + k*SZ1*SZ2] = Mod.O.C_ext;
            }  
        }
    }

    /*doing a cross section then repeating over the length*/
    for(i=0;i<SZ1;i++)
    {
        /*placing hepatocyte*/
        /*filling first and last line*/
        for(j=0;j<SZ2;j++)
        {  
            for(k=0;k<10;k++)
            {
               place_hepatocyte(i,j,k);  
            }

            for(k=SZ3-10;k<SZ3;k++)
            {
               place_hepatocyte(i,j,k);  
            }
        }


        for(k=10;k<(SZ3-10);k++)
        {
            
            for(j=0;j<10;j++)
            {
                place_hepatocyte(i,j,k);
            }
            
            for(j=SZ2-10;j<SZ2;j++)
            {
               place_hepatocyte(i,j,k);
            }
        }
        
        /*the space of Disse is filled with matrix protein*/
        /*just above the hepatocyte layer*/
        for(j=10;j<SZ2-10;j++)
        {
            place_Disse(i,j,10);
            place_Disse(i,j,SZ3-11);
        }

        for(k=11;k<SZ3-11;k++)
        {
            place_Disse(i,10,k);
            place_Disse(i,SZ2-11,k);
        }

        /*the lsec layer is placed*/
        /*just above the hepatocyte layer*/
        for(j=11;j<SZ2-11;j++)
        {
            place_lsec(i,j,11);
            place_lsec(i,j,SZ3-12);
        }

        for(k=12;k<SZ3-12;k++)
        {
            place_lsec(i,11,k);
            place_lsec(i,SZ3-12,k);
        }
    }

    /*place stellate cells in the space of Disse*/
    place_stellate(5,10,SZ3/2);
    place_stellate(25,SZ2/2,10);
    place_stellate(45,SZ2-11,SZ3/2);
    place_stellate(65,SZ2/2,SZ3-11);
    place_stellate(85,10,SZ3/2);


    place_Kupffer(10,SZ2/2,SZ3/2);
    place_Kupffer(30,SZ2/2,SZ3/2);
    place_Kupffer(50,SZ2/2,SZ3/2);
    place_Kupffer(70,SZ2/2,SZ3/2);
    place_Kupffer(90,SZ2/2,SZ3/2);

}

