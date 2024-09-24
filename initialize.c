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
    Mod.S.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.S.DC_tiss;
    Mod.S.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kS;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    /*updating the cell array*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 1; /*hepatocyte cell type=1*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_cons = Mod.kS;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_diff = Mod.S.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_diff = 0.0;
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
    Mod.S.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.S.DC_mat;
    Mod.S.kC[row + col*SZ1 +stck*SZ1*SZ2] = 0.0;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_mat;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = 0.0;
}

void place_lsec(int row, int col, int stck)
{
    /*increasing cell number*/
    Mod.M_Tissue.N_Cell++;

    /*placing things on Grids*/
    Mod.M_Tissue.Grid[row + col*SZ1 +stck*SZ1*SZ2] =  Mod.M_Tissue.N_Cell;
    Mod.S.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.S.DC_endo;
    Mod.S.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kS;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_endo;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 2; /*lse cell type=2*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_cons = Mod.kS;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_diff = Mod.S.DC_endo;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_diff = 0.0;
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
    Mod.S.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.S.DC_tiss;
    Mod.S.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kS;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 3; /*stellate cell type=3*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_cons = Mod.kS;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_diff = Mod.S.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = Mod.kO;
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
    Mod.S.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.S.DC_tiss;
    Mod.S.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kS;
    Mod.O.DCm[row + col*SZ1 +stck*SZ1*SZ2] = Mod.O.DC_tiss;
    Mod.O.kC[row + col*SZ1 +stck*SZ1*SZ2] = Mod.kO;

    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_index = Mod.M_Tissue.N_Cell;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cell_type= 4; /*Kupffer cell type=4*/
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].x_pos= row; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].y_pos= col; 
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].z_pos= stck;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_cons = Mod.kS;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].S_diff = Mod.S.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dS_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_cons = Mod.kO;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_cons = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].O_diff = Mod.O.DC_tiss;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].dO_diff = 0.0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].parent_idx = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Chg_timer = 0;
    Mod.M_Cell[Mod.M_Tissue.N_Cell-1].Cycle_dur = 1500;
}

void initialize(float kS, float kO)
{
    printf("initializing... \n");

    Mod.kO = kO;
	Mod.kS = kS;

    int i=0,j=0,k=0,l=0;

    /*setting constants*/
    Mod.M_Tissue.N_Cell =0;
    Mod.S.DC_med = 40000.0;
    Mod.S.DC_mat = 30000.0;
    Mod.S.DC_tiss = 7000.0;
    Mod.S.DC_endo = 3000.0;
    Mod.O.DC_med = 200000.0;
    Mod.O.DC_mat = 180000.0;
    Mod.O.DC_tiss = 120000.0;
    Mod.O.DC_endo = 120000.0;


    /*allocating memory for the nutrient grids*/
    Mod.S.C = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.S.DCm = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.S.kC = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.S.fl = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.S.dx = 2.0;

    Mod.O.C = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.DCm = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.kC = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.fl = malloc(SZ1*SZ2*SZ3*sizeof(double));
    Mod.O.dx = 2.0;

    /*Grid to follow cells*/
    Mod.M_Tissue.Grid = malloc(SZ1*SZ2*SZ3*sizeof(int));
	Mod.M_Tissue.LD = malloc(SZ1*SZ2*SZ3*sizeof(int));
	Mod.M_Tissue.state_mat = malloc(SZ1*SZ2*SZ3*sizeof(unsigned char));
	Mod.M_Cell = malloc(100000*sizeof(Mod.M_Cell));


    for(k=0;k<SZ3;k++)
    {
        for(j=0;j<SZ2;j++)
        {
              for(i=0;i<SZ1;i++)
            {
                Mod.S.DCm[i + j*SZ1 + k*SZ1*SZ2] = Mod.S.DC_med;
                Mod.O.DCm[i + j*SZ1 + k*SZ1*SZ2] = Mod.O.DC_med;
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
            place_hepatocyte(i,j,0);
            place_hepatocyte(i,j,SZ3-1);
        }

        for(k=1;k<SZ3-1;k++)
        {
            place_hepatocyte(i,0,k);
            place_hepatocyte(i,SZ2-1,k);
        }
        
        /*the space of Disse is filled with matrix protein*/
        /*just above the hepatocyte layer*/
        for(j=1;j<SZ2-2;j++)
        {
            place_Disse(i,j,1);
            place_Disse(i,j,SZ3-2);
        }

        for(k=2;k<SZ3-2;k++)
        {
            place_Disse(i,1,k);
            place_Disse(i,SZ2-2,k);
        }

        /*the lsec layer is placed*/
        /*just above the hepatocyte layer*/
        for(j=2;j<SZ2-2;j++)
        {
            place_lsec(i,j,2);
            place_lsec(i,j,SZ3-3);
        }

        for(k=3;k<SZ3-3;k++)
        {
            place_lsec(i,2,k);
            place_lsec(i,SZ3-3,k);
        }
    }

    /*place stellate cells in the space of Disse*/
    place_stellate(5,1,7);
    place_stellate(25,7,1);
    place_stellate(45,SZ2-2,7);
    place_stellate(65,7,SZ3-2);
    place_stellate(85,1,7);


    place_Kupffer(10,7,7);
    place_Kupffer(30,7,7);
    place_Kupffer(50,7,7);
    place_Kupffer(70,7,7);
    place_Kupffer(90,7,7);

}

