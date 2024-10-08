#define SZ1 100
#define SZ2 35
#define SZ3 35

    typedef struct Nutrs
    {
		/*nutrient*/
		double *C;

		/*Diffusion*/
		double *DCm;
		double DC_med;
		double DC_tiss;
		double DC_mat;
		double DC_endo;
		
		/*ext conc.*/
		double C_ext;
		
		/*consommation*/
		double *kC;
		double cC; /*parameter for substrate Hill function*/

		/*flow*/
		double *fl_x;
	
		/* model parameters*/
		double dx;
		double dt;

    }Nutr;

	typedef struct Tissues
	{
		/* Cell & State Grids*/
		int *Grid;
		int *LD;
		unsigned char *state_mat;

        double G_hypo;
		double O_norm;


        /*other variables*/
		int N_Cell; /*all cells (dead or alive) in the model*/
		int N_Dead;
		int N_Live;
		int reac_time; /*time for consumption changes*/
		int avr_cell_cyc; /*cell cycle reaction*/
		int n_pts;
		int mean_rad_sq;
		int ctr_x;
		int ctr_y;
		int ctr_z;
		int rad_pellet;/*pellet radius for chip model*/
    }Tissue;

    typedef struct Cells
	{
		int Cell_index;
		int Timer;
		int Cycle_dur;
		int Chg_timer;
		float G_cons;
		float dG_cons;
		float G_diff;
		float dG_diff;
		float O_cons;
		float dO_cons;
		float O_diff;
		float dO_diff;
		int x_pos;
		int y_pos;
		int z_pos;
		int parent_idx;
		int Cell_type;

	}Cell;


    typedef struct Models
	{
		Cell *M_Cell;
		Tissue M_Tissue;
		Nutr G;
        Nutr O;
		Nutr T;
		int elapsed_mins;
		int n_pts;
		double kG;/*model consumptions needed to restart calculations*/
		double kO; 
        
	}Model;