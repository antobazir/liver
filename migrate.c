#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structdef.h"

extern Model Mod;

void shift(int cell_idx, int row, int col, int stck, int new_row, int new_col, int new_stck)
{
    				/*duplication*/
				Mod.M_Tissue.Grid[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = cell_idx+1;
				Mod.M_Tissue.LD[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = 1;
				Mod.S.DCm[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = Mod.S.DCm[row + col * SZ1 + stck * SZ1 * SZ2];
				Mod.O.DCm[new_row + new_col*SZ1 + new_stck*SZ1*SZ2]= Mod.O.DCm[row + col * SZ1 + stck * SZ1 * SZ2];
				Mod.S.kC[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = Mod.S.kC[row + col * SZ1 + stck * SZ1 * SZ2];
				Mod.O.kC[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = Mod.O.kC[row + col * SZ1 + stck * SZ1 * SZ2];
				Mod.M_Tissue.state_mat[new_row + new_col*SZ1 + new_stck*SZ1*SZ2] = Mod.M_Tissue.state_mat[row + col * SZ1 + stck * SZ1 * SZ2];
				Mod.M_Cell[cell_idx].x_pos = new_row;
				Mod.M_Cell[cell_idx].y_pos = new_col;
				Mod.M_Cell[cell_idx].z_pos = new_stck;


				/*removal*/
				Mod.M_Tissue.Grid[row + col * SZ1 + stck * SZ1 * SZ2] = 0;
				Mod.M_Tissue.LD[row + col * SZ1 + stck * SZ1 * SZ2] = 0;
				if(Mod.M_Cell[cell_idx].Cell_type==3)
				{	
					Mod.S.DCm[row + col * SZ1 + stck * SZ1 * SZ2] = Mod.S.DC_mat;
					Mod.O.DCm[row + col * SZ1 + stck * SZ1 * SZ2] = Mod.O.DC_mat;
				}
				
				if(Mod.M_Cell[cell_idx].Cell_type==4)
				{	
					Mod.S.DCm[row + col * SZ1 + stck * SZ1 * SZ2] = Mod.S.DC_med;
					Mod.O.DCm[row + col * SZ1 + stck * SZ1 * SZ2] = Mod.O.DC_med;
				}
				Mod.S.kC[row + col * SZ1 + stck * SZ1 * SZ2] = 0.0;
				Mod.O.kC[row + col * SZ1 + stck * SZ1 * SZ2] = 0.0;
				Mod.M_Tissue.state_mat[row + col * SZ1 + stck * SZ1 * SZ2] = 0;
}

void migrate(int cell_index)
{
    /*stellate cells migrate if they are activated and stays in the space of Disse*/

    int i, j, k, l, n, m;
	int row_perim[10000];
	int col_perim[10000];
	int stck_perim[10000];
	int n_row[26];
	int n_col[26];
	int n_stck[26];
	int nf_row[26];
	int nf_col[26];
	int nf_stck[26];
	int row_shift[26];  
	int col_shift[26];
	int stck_shift[26];
	int idx_ngh[26];
	int row_ctrd;
	int col_ctrd;
	int stck_ctrd;
	int lgth_perim = 0;
	int rad_sq[10000];
	int rad_sq_n[26];
	int min_rad;
	int n_closer, n_closest;
	int closer[26], closest[26];
	int chosen_site;
	int n_free;
	int row, col, stck;
	int ngh[3][3][3];
	double allowed_env;
	/*FILE *f;*/

	m=0;
	n_free=0;
	n=0;

	for (i = -1; i < 2; i++)
	{
		for (j = -1; j < 2; j++)
		{
			for (m = -1; m < 2; m++)
			{
				/*fais tous les voisins en évitant le point central*/
				if (i != 0 || j != 0 || m != 0)
				{
					row_shift[l] = i;
					col_shift[l] = j;
					stck_shift[l] = m;
					n++;
				}
			}
		}
	}

	m=0;
	n=0;
	srand(cell_index);

    row = Mod.M_Cell[cell_index].x_pos;
	col = Mod.M_Cell[cell_index].y_pos;
	stck = Mod.M_Cell[cell_index].z_pos;

	/*if the cell moving is not on the borders...*/
	if((row>0)&&(row<SZ1-1))
	{
		if (Mod.M_Cell[cell_index].Cycle_dur > 1.01)
		{
			/*stellate cells can only mode in matrix*/
			if(Mod.M_Cell[cell_index].Cell_type==3)
			{
				allowed_env = Mod.S.DC_mat;
			}
		
			/*Kupffer cells can only mode in blood/medium*/
			if(Mod.M_Cell[cell_index].Cell_type==4)
			{
				allowed_env = Mod.S.DC_med;
			}

			for(l = -1; l<2 ; l++)
			{ 
				for (j = -1; j < 2; j++)
				{
					for (i = -1; i < 2; i++)
					{
						ngh[i + 1][j + 1][l + 1] = Mod.S.DCm[row + i + (col + j)*SZ1 + (stck + l)*SZ1*SZ2];
						/*printf("current: %d \n",ngh[i+1][j+1]);*/
						n_row[m] = row + i;
						n_col[m] = col + j;
						n_stck[m] = stck + l;
						/*if the neighbor is matrix, it is open*/
						if (ngh[i + 1][j + 1][l + 1] == allowed_env)
						{
							if (i != 0 || j != 0||l != 0 )
							{
								nf_row[m] = row + i;
								nf_col[m] = col + j;
								nf_stck[m] = stck + l;
								idx_ngh[m] = n + 1;
								n_free = n_free + 1;
								m = m + 1;
							}

						}
						n =n + 1;
					}
				}
			}
		

			if (n_free == 0)
			{
				/*can't move if boxed*/
				return;
			}

			/* one free neighbor*/
			if (n_free == 1)
			{
				/*movement =duplication + removal*/

				/*duplication*/
				
				shift(cell_index,row,col,stck,nf_row[0],nf_col[0],nf_stck[0]);
			}

			if (n_free > 1)
			{
				chosen_site = (rand() %
					((n_free - 1) - 0 + 1)) + 0;

				/*printf("chosen_site: %d\n",chosen_site);*/
							/*movement =duplication + removal*/
				shift(cell_index,row,col,stck,nf_row[chosen_site],nf_col[chosen_site],nf_stck[chosen_site]);
			}
		}
	}

	else
	{
		if((row==0))
		{
			/*shift back to "safe" spot or do nothing and wait for next round*/
			if(Mod.M_Tissue.Grid[row+1 + col*SZ1 + stck*SZ1*SZ2]==0)
			{
				shift(cell_index,row,col,stck,row+1,col,stck);
			}
		}

		if((row==SZ1-1))
		{
			/*shift back to "safe" spot or do nothing and wait for next round*/
			if(Mod.M_Tissue.Grid[row-1 + col*SZ1 + stck*SZ1*SZ2]==0)
			{
				shift(cell_index,row,col,stck,row-1,col,stck);
			}
		}
	}
    return;
}
