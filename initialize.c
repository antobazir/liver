#include<stdio.h>
#include"structdef.h"/*the header containing model types*/

extern Model Mod;/*line that link this invokes the model instance called in liver.c*/

void initialize(float kS, float kO,float behav_code)
{
    printf("initializing... \n");

    Mod.kO = kO;
	Mod.kS = kS;

    int i=0,j=0,k=0,l=0;

    /*allocating memory for the nutrient grids*/
    Mod.S.C = malloc(SZ1*SZ2*SZ3,sizeof(double));
    Mod.S.DCm = malloc(SZ1*SZ2*SZ3,sizeof(double));
    Mod.S.kC = malloc(SZ1*SZ2*SZ3,sizeof(double));
    Mod.O.C = malloc(SZ1*SZ2*SZ3,sizeof(double));
    Mod.O.DCm = malloc(SZ1*SZ2*SZ3,sizeof(double));
    Mod.O.kC = malloc(SZ1*SZ2*SZ3,sizeof(double));

    Mod.M_Tissue.Grid = malloc(SZ1*SZ2*SZ3,sizeof(int));
	Mod.M_Tissue.LD = malloc(SZ1*SZ2*SZ3,sizeof(int));
	Mod.M_Tissue.state_mat = malloc(SZ1*SZ2*SZ3,sizeof(unsigned char));
	Mod.M_Cell = malloc(100000*sizeof(Mod.M_Cell));
}