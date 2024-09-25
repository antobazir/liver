#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structdef.h"

/*This modules runs the diffusion-advection-reaction model for all molelcules*/
/*The flow is only in the x-direction*/

extern Model Mod;



void diff_adv_reac()
{

  omp_set_num_threads(4);
	int i,j,k,l,m,h,ntS;
	ntS = (int)(1.0/Mod.S.dt);
	int rad_sq;

  for (h=0; h<ntS; h++)
	{
    /*openmp directive to disitribute work across cores*/
		#pragma omp parallel for collapse(3)
		for (m=1;m<SZ3-1;m++)
		{
			for (j=1; j<SZ2-1; j++)
			{	
				for(i=1; i<SZ1-1; i++)
				{
					Mod.S.C[i + j*SZ1 + m*SZ1*SZ2] = Mod.S.C[i + j*SZ1 + m*SZ1*SZ2] + 

          /*upwind x-direction flow term*/
					Mod.S.fl_x[i + j*SZ1 + m*SZ1*SZ2]*Mod.S.dt/Mod.S.dx*(Mod.S.C[(i) + j*SZ1 + m*SZ1*SZ2]-Mod.S.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+

          /*3D central difference explicit diffusion terms w/ variable diffusivity*/
          Mod.S.dt/(Mod.S.dx*Mod.S.dx)*(
						(Mod.S.DCm[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.S.DCm[(i-1) + j*SZ1 + m*SZ1*SZ2])*(Mod.S.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.S.C[(i-1) + j*SZ1 + m*SZ1*SZ2])/4+
						(Mod.S.DCm[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.S.DCm[i+ (j-1)*SZ1 + m*SZ1*SZ2])*(Mod.S.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.S.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])/4+
						(Mod.S.DCm[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.S.DCm[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])*(Mod.S.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.S.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])/4+
						Mod.S.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.S.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-2*Mod.S.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.S.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+
						Mod.S.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.S.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-2*Mod.S.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.S.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])+
						Mod.S.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.S.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-2*Mod.S.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.S.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])
						)

            /*consumption term in form of a Hill function*/
						-Mod.S.dt*Mod.S.kC[i + j*SZ1 + m*SZ1*SZ2]*(Mod.S.C[i + j*SZ1 + m*SZ1*SZ2]/(Mod.S.C[i + j*SZ1 + m*SZ1*SZ2]+0.05));

          /*little fail-safe that shouldn't be needed*/
					if(Mod.S.C[i + j*SZ1 + m*SZ1*SZ2]<0.0)
						{
							Mod.S.C[i + j*SZ1 + m*SZ1*SZ2]==0.0;
						} 
					
				}	

			}

		}


			/*Boundary Conditions*/

      /*x-axis*/
      #pragma omp parallel for collapse(2)
			for (m=0;m<SZ3;m++)
			{
				for (l=0; l<SZ2; l++)
				{
          /*dirichlet constant value on the inlet/portal field side*/
					Mod.S.C[0 + l*SZ1 + m*SZ1*SZ2] = Mod.S.C_ext ;

          /*outgoing flux on the outlet/central vein side*/
					Mod.S.C[SZ1-1 + l*SZ1 + m*SZ1*SZ2]  = Mod.S.dt/(Mod.S.dx*Mod.S.dx)*(Mod.S.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.S.C[(SZ1-1) + j*SZ1 + m*SZ1*SZ2]-2*Mod.S.C[(SZ1-2) + j*SZ1 + m*SZ1*SZ2] +Mod.S.C[(SZ1-3) + j*SZ1 + m*SZ1*SZ2]));
				}
				}
			}

      /*y-axis*/
			#pragma omp parallel for collapse(2)
			for (m=0;m<SZ3;m++)
			{
				for(k=0; k<SZ1; k++)
				{		
          /*zero flux conditions on both sides*/
					Mod.S.C[k + 0*SZ1 + m*SZ1*SZ2] = Mod.S.C[k + 1*SZ1 + m*SZ1*SZ2] ;
					Mod.S.C[k + (SZ2-1)*SZ1 + m*SZ1*SZ2] = Mod.S.C[k + (SZ2-2)*SZ1 + m*SZ1*SZ2];
			}

      /*z-axis*/
			#pragma omp parallel for collapse(2)
			for (l=0; l<SZ2; l++)
			{
				for(k=0; k<SZ1; k++)
				{	
          /*zero flux conditions on both sides*/
					Mod.S.C[k + l*SZ1 + 0*SZ1*SZ2] = Mod.S.C[k + l*SZ1 + 1*SZ1*SZ2] ;
					Mod.S.C[k + l*SZ1 + (SZ3-1)*SZ1*SZ2]= Mod.S.C[k + l*SZ1 + (SZ3-2)*SZ1*SZ2];
				}
      }	
		
	}

      return;
}