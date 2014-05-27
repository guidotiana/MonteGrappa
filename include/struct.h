/*
 * struct.h
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

struct vector { double x,y,z; };

struct angles { double ang;
                double dih;
                double r;
              };


struct s_polymer
{
	int nback;			// number of backbone atoms
	char title[80];
	struct s_back *back;
	struct vector **vback;		// lookback table (given an atom, gives the address of the cartesian position; set by ReadPolymer)
	struct s_tables *tables;	// tables of sin,cos,etc. ( needs InitTables before use)
	double etot;			// total energy
	double t;			// its temperature (necessary for op)
      double **A,**G,**L,**Y;      //local move
      double *g_ang,*d_ang;        //loval move



	#ifdef OPTIMIZEPOT
	struct s_optimizepot *op;
	#endif
};

struct s_back
{
	int ia;
	int iapdb;
	int itype;
	char aa[5];
	char type[5];
	int iaa;
	struct vector pos;		// cartesian position
	struct angles sph;		// position in spherical coordinates based on the 3 preceding backbone atoms
	int nside;			// number of atoms in its sidechain (0=none)
	int move;			// 0 means that it should't be moved
	int irot;			// current rotamer id
	int nrot;			// total number of rotamers

	struct s_side *side;		// sidechain atoms

	int ncontacts;			// number of backbone atoms whose atom (including sidechain) are in contact with this
	int *contacts;			// list of atoms in contact
	int *contacts_p;		// ... which chain it belongs to
	double *e;			// interaction energy of this with the others residues
	double e_ang;			// angular energy
	double e_dih;			// dihedral energy
	int nshell;			// number of backbone atoms in the shell
	int *shell;			// atoms in the shell
	int *shell_p;			// ... which chain it belongs to

	double d2_next;			// square distance to the next backbone
	double a_next;			// angle
	double d_next;			// dihedral
};

struct s_side
{
	int ia;
	int iapdb;
	int itype;
	char type[5];
	struct vector pos;		// cartesian position

	struct s_rotamers *rot;		// rotamers
};

struct s_rotamers
{
	int b1;						// the atoms which act as reference state
	int b2;
	int b3;
	struct angles ang;				// the spherical coordinate of the sidechain atom
};

struct s_potential
{
	double **e;					// energy of the square well
	double **r_2;					// square width of the square well
	double **r0_2;					// square width of hard core repulsion

	double *e_ang;					// energy constant of angular potential
	double *ang0;					// equilibrium angle

	int dih_periodic;				// 1=activate periodic dihedrals
	double *e_dih1;					// energy constant of dihedral potential with multeplicity 1
	double *dih01;					// phase shift of diherdrals with multeplicity 1
	double *e_dih3;					// energy constant of dihedral potential with multeplicity 3
	double *dih03;					// phase shift of diherdrals with multeplicity 3

	int dih_tabled;					// 1=activate tabled dihedrals
	int *dih_which;					// 0=phi, 1=psi
	double *dih_pa;					// weight of alpha-helix in tabled dihedral
	double *dih_pb;					// weight of beta-sheet in tabled dihedrals
	double *dih_f_psi_a;				// terms of the potential
	double *dih_f_phi_a;
	double *dih_f_psi_b;
	double *dih_f_phi_b;

	double g_r0hard;				// global hardcore repulsion
	double g_ehomo;					// global homopolymer e
	double g_rhomo;
	double g_anglek;
	double g_angle0;
	double g_dihk1;
	double g_dih01;
	double g_dihk3;
	double g_dih03;
	int g_imin;					// imin
	double g_dihbin;				// bin of tabled dihedral potential
	double g_dihe;					// energy for tabled dihedral potential
	char boxtype;					// n=none, c=cubic, s=spherical
	double boxsize;

	double lbox;

	double kr2_splice;				// splice energy well into two parts, the distance being kr * r
	double ke_splice;				// and the depth is ke * e
	int splice;						// =1 to activate

	int *hc_type;					// hardcore parameters
	double *hc_r0;
	int hc_number;

	int *hb;					// hydrogen bond: 0=none, 1=donor, 2=acceptor
	int *hb_iam;					// ia of the atom preceding donor/acceptor

	double **sigma;				//energy dihedrals with minimum in Ramachandran dihedrals
	int **dih0;
	int dih_ram;
	double e_dihram;
	double **ab_propensity;
};

struct s_mc_parms
{
	int npol;		// number of chains
	int nstep;
	long seed;
	double dw_flip;		// maximum angle allowed in a flip
	double dw_pivot;	// maximum angle allowed in a flip
	double dw_mpivot;	// maximum angle allowed in a flip
	double dw_lpivot;	// maximum angle allowed in a flip
	double dw_mflip;	// maximum angle allowed in a multiple flip
	char fntrj[50];		// name of trajectory file
	char fne[50];		// name of energy file
	char flastp[50];	// name of last conformation;
	char fnproc[50];	// name of the log file
	int nprinttrj;		// when to print trajectory file
	int nprintlog;		// when to print log
	int nprinte;		// when to print energy
	FILE *flog;
	int shell;		// 1=activate shell
	int nshell;		// update shell every %d step
	double r2shell;		// square distance to be in the shell of a backbone atom
	int ntemp;		// number of temperatures (replicas)	
	#ifdef ACTIVE_MPI
	double T[NTEMPMAX];	// temperatures
	#else
	double T;
	#endif
	int randdw;		// 1=dw from flat distribution, 2=dw from gaussian distribution
	int debug;
	int movetype[NMOVES];	//activate(0/1) move: 0=flip, 1=pivot, 2=multipivot, 3=sidechain, 4=loose multipivot
	int nmul_mpivot;	// how many backbone atoms to move in multiple pivot move
	int nmul_lpivot;	// how many backbone atoms to move in loose pivot move
	int nmul_mflip;		// how many backbone atoms to move at most in multiple flip move
	int nosidechains;	//0=there are, 1=no sidechains
	int noangpot;		//0=there is a potential on angles, 1; there is not
	int nodihpot;		//0=there is a potential on dihedrals, 1; there is not
	int nrun;		// number of repetitions of the MC
	int always_restart;	//1=in different irun, start always from input structure
	int record_native; //1=records the input (native( structure as first snapshot
	int acc;		// # of accepted moves
	int mov;		// # of attempted moves
	int disentangle;	//0=if two atoms are closer than r_hardcore, then do not calculate the energy
	int stempering;		//1=do stempering
	double dx_com;		//displacemente oc center of mass
	double dx_clm;		//displacement of cluster
	double r_cloose;	// constrains in the bond distance,
	double a_cloose;	// angles
	double d_cloose;	// and dihedral, relative to the initial position (-1 to disable)
	int hb;			// activate hydrogen bonds
	#ifdef OPTIMIZEPOT
	 struct s_optimizepot_input *op_input;
	 char fnop[50];		// name of file of experimental restrains
	 char op_minim[50];	// optimize procedure: none, sample, steep
	 int op_itermax;	// # of iteration steps
	 double op_step;
	 double op_T;		// temperature corresponding to the experimental data (can be different from T)
	 int op_deltat;		// = nstep/op_nframes
	 double op_stop;	// stop condition
	 int op_print;		// print log every * step
	 double op_emin;	// matrix elements cannot go below this limit
	 double op_emax;	// matrix elements cannot go below this limit
	 int op_wait;		// discard the first steps
	 double op_r;		// default well width;
	 double op_r0;		// default well hardcore;
	#endif
	int anneal;		// 1=anneal,0=not 
	int anneal_often;	// every * istep 
	int anneal_step;	// carry out * step at higher temperature (* decreases to zero)
	double anneal_t;	// higher temperature
	int anneal_recov;	// do nothing for * steps after returning to actual temperature
	

	int nstep_exchange;	
	int nmul_local;
};

struct s_tables
{
	double *fast_sqrt;
	double *fast_sin;
	double *fast_cos;
	double *fast_expp;		// exp of positive numbers
	double *fast_expm;		// exp of negative numbers
	double *fast_acos;		// arccos

};
