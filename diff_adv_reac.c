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
  omp_set_num_threads(2);
	int i,j,k,l,m,h,ntG,ntO;
	ntG = (int)(60.0/Mod.O.dt);
	int rad_sq;

  for (h=0; h<ntG; h++)
	{
    /*openmp directive to disitribute work across cores*/
		#pragma omp parallel for collapse(3)
		for (m=1;m<SZ3-1;m++)
		{
			for (j=1; j<SZ2-1; j++)
			{	
				for(i=1; i<SZ1-1; i++)
				{
					Mod.G.C[i + j*SZ1 + m*SZ1*SZ2] = Mod.G.C[i + j*SZ1 + m*SZ1*SZ2] -

          /*upwind x-direction flow term*/
					Mod.G.fl_x[j + m*SZ2]*Mod.G.dt/Mod.G.dx*(Mod.G.C[(i) + j*SZ1 + m*SZ1*SZ2]-Mod.G.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+

          /*3D central difference explicit diffusion terms w/ variable diffusivity*/
          Mod.G.dt/(Mod.G.dx*Mod.G.dx)*(
						(Mod.G.DCm[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.G.DCm[(i-1) + j*SZ1 + m*SZ1*SZ2])*(Mod.G.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.G.C[(i-1) + j*SZ1 + m*SZ1*SZ2])/4+
						(Mod.G.DCm[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.G.DCm[i+ (j-1)*SZ1 + m*SZ1*SZ2])*(Mod.G.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.G.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])/4+
						(Mod.G.DCm[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.G.DCm[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])*(Mod.G.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.G.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])/4+
						Mod.G.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.G.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-2*Mod.G.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.G.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+
						Mod.G.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.G.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-2*Mod.G.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.G.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])+
						Mod.G.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.G.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-2*Mod.G.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.G.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])
						)

            /*consumption term in form of a Hill function*/
						-Mod.G.dt*Mod.G.kC[i + j*SZ1 + m*SZ1*SZ2]*(Mod.G.C[i + j*SZ1 + m*SZ1*SZ2]/(Mod.G.C[i + j*SZ1 + m*SZ1*SZ2]+0.05));

          /*little fail-safe that shouldn't be needed*/
					if(Mod.G.C[i + j*SZ1 + m*SZ1*SZ2]<0.0)
					{
							Mod.G.C[i + j*SZ1 + m*SZ1*SZ2]==0.0;
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

        if(Mod.G.fl_x[l + m*SZ2]!=0)
        {
          /*dirichlet constant value on the inlet/portal field side*/
          Mod.G.C[0 + l*SZ1 + m*SZ1*SZ2] = Mod.G.C_ext ;

          /*outgoing flux on the outlet/central vein side*/
          Mod.G.C[SZ1-1 + l*SZ1 + m*SZ1*SZ2]  = 2*Mod.G.C[SZ1-2 + l*SZ1 + m*SZ1*SZ2] - Mod.G.C[SZ1-3 + l*SZ1 + m*SZ1*SZ2];
        }

        else
        {
          /*inside the hepatocyte zero flux on borders*/
          Mod.G.C[0 + l*SZ1 + m*SZ1*SZ2] = Mod.G.C[1 + l*SZ1 + m*SZ1*SZ2];
          Mod.G.C[SZ1-1 + l*SZ1 + m*SZ1*SZ2] = Mod.G.C[SZ1-2 + l*SZ1 + m*SZ1*SZ2] ;
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
        Mod.G.C[k + 0*SZ1 + m*SZ1*SZ2] = Mod.G.C[k + 1*SZ1 + m*SZ1*SZ2] ;
        Mod.G.C[k + (SZ2-1)*SZ1 + m*SZ1*SZ2] = Mod.G.C[k + (SZ2-2)*SZ1 + m*SZ1*SZ2];
      }
    }

    /*z-axis*/
    #pragma omp parallel for collapse(2)
    for (l=0; l<SZ2; l++)
    {
      for(k=0; k<SZ1; k++)
      {	
        /*zero flux conditions on both sides*/
        Mod.G.C[k + l*SZ1 + 0*SZ1*SZ2] = Mod.G.C[k + l*SZ1 + 1*SZ1*SZ2] ;
        Mod.G.C[k + l*SZ1 + (SZ3-1)*SZ1*SZ2]= Mod.G.C[k + l*SZ1 + (SZ3-2)*SZ1*SZ2];
      }
    }	
		/*printf("%f\n",Mod.G.C[SZ1/2 + SZ2/2*SZ1 + SZ3/2*SZ1*SZ2]);*/
	}

  /*Oxygen solving*/
  for (h=0; h<ntO; h++)
	{
    /*openmp directive to disitribute work across cores*/
		#pragma omp parallel for collapse(3)
		for (m=1;m<SZ3-1;m++)
		{
			for (j=1; j<SZ2-1; j++)
			{	
				for(i=1; i<SZ1-1; i++)
				{
					Mod.O.C[i + j*SZ1 + m*SZ1*SZ2] = Mod.O.C[i + j*SZ1 + m*SZ1*SZ2] -

          /*upwind x-direction flow term*/
					Mod.G.fl_x[j + m*SZ2]*Mod.O.dt/Mod.O.dx*(Mod.O.C[(i) + j*SZ1 + m*SZ1*SZ2]-Mod.O.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+

          /*3D central difference explicit diffusion terms w/ variable diffusivity*/
          Mod.O.dt/(Mod.O.dx*Mod.O.dx)*(
						(Mod.O.DCm[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.O.DCm[(i-1) + j*SZ1 + m*SZ1*SZ2])*(Mod.O.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-Mod.O.C[(i-1) + j*SZ1 + m*SZ1*SZ2])/4+
						(Mod.O.DCm[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.O.DCm[i+ (j-1)*SZ1 + m*SZ1*SZ2])*(Mod.O.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-Mod.O.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])/4+
						(Mod.O.DCm[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.O.DCm[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])*(Mod.O.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-Mod.O.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])/4+
						Mod.O.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.O.C[(i+1) + j*SZ1 + m*SZ1*SZ2]-2*Mod.O.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.O.C[(i-1) + j*SZ1 + m*SZ1*SZ2])+
						Mod.O.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.O.C[i+ (j+1)*SZ1 + m*SZ1*SZ2]-2*Mod.O.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.O.C[i+ (j-1)*SZ1 + m*SZ1*SZ2])+
						Mod.O.DCm[i + j*SZ1 + m*SZ1*SZ2]*(Mod.O.C[i+ (j)*SZ1 + (m+1)*SZ1*SZ2]-2*Mod.O.C[i + j*SZ1 + m*SZ1*SZ2] +Mod.O.C[i+ (j)*SZ1 + (m-1)*SZ1*SZ2])
						)

            /*consumption term in form of a Hill function*/
						-Mod.O.dt*Mod.O.kC[i + j*SZ1 + m*SZ1*SZ2]*(Mod.O.C[i + j*SZ1 + m*SZ1*SZ2]/(Mod.O.C[i + j*SZ1 + m*SZ1*SZ2]+0.05));

          /*little fail-safe that shouldn't be needed*/
					if(Mod.O.C[i + j*SZ1 + m*SZ1*SZ2]<0.0)
					{
							Mod.O.C[i + j*SZ1 + m*SZ1*SZ2]==0.0;
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
        if(Mod.G.fl_x[l + m*SZ2]!=0)
        {
          /*dirichlet constant value on the inlet/portal field side*/
          Mod.O.C[0 + l*SZ1 + m*SZ1*SZ2] = Mod.O.C_ext ;

          /*outgoing flux on the outlet/central vein side*/
          Mod.O.C[SZ1-1 + l*SZ1 + m*SZ1*SZ2]  = 2*Mod.O.C[SZ1-2 + l*SZ1 + m*SZ1*SZ2] - Mod.O.C[SZ1-3 + l*SZ1 + m*SZ1*SZ2];
        }

        else
        {
          /*inside the hepatocyte zero flux on borders*/
          Mod.O.C[0 + l*SZ1 + m*SZ1*SZ2] = Mod.O.C[1 + l*SZ1 + m*SZ1*SZ2];
          Mod.O.C[SZ1-1 + l*SZ1 + m*SZ1*SZ2] = Mod.O.C[SZ1-2 + l*SZ1 + m*SZ1*SZ2] ;
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
        Mod.O.C[k + 0*SZ1 + m*SZ1*SZ2] = Mod.O.C[k + 1*SZ1 + m*SZ1*SZ2] ;
        Mod.O.C[k + (SZ2-1)*SZ1 + m*SZ1*SZ2] = Mod.O.C[k + (SZ2-2)*SZ1 + m*SZ1*SZ2];
      }
    }

    /*z-axis*/
    #pragma omp parallel for collapse(2)
    for (l=0; l<SZ2; l++)
    {
      for(k=0; k<SZ1; k++)
      {	
        /*zero flux conditions on both sides*/
        Mod.O.C[k + l*SZ1 + 0*SZ1*SZ2] = Mod.O.C[k + l*SZ1 + 1*SZ1*SZ2] ;
        Mod.O.C[k + l*SZ1 + (SZ3-1)*SZ1*SZ2]= Mod.O.C[k + l*SZ1 + (SZ3-2)*SZ1*SZ2];
      }
    }	
		/*printf("%f\n",Mod.G.C[SZ1/2 + SZ2/2*SZ1 + SZ3/2*SZ1*SZ2]);*/
	}
  
  return;
}