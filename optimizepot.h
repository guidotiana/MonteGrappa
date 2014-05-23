/*
 * optimizepot.h
 *
 *  Created on: Jan 31, 2011
 *      Author: guido
 */



#define OP_NCONTMAX	500	// maximum number of contacts
#define OP_NRESMAX	100	// maximum number of restrains
#define REALLOCATE		// if the number of contacts exceed OP_NCONTMAX, reallocate the structure
//#define OP_DEBUG

struct s_optimizepot *InitializeOptimizePot(struct s_mc_parms *parms, int ntypes, struct s_potential *u, FILE *fproc, int iproc, FILE *fp, int *nrestr);
struct s_optimizepot *AlloOptimizePot(struct s_mc_parms *parms, int ntypes, int nrestr, FILE *fproc);
struct s_optimizepot_input *ReadOPRestrains(struct s_mc_parms *parms);
void OP_GetRestrain(int it, struct s_polymer *p, int ipol, struct s_optimizepot_input *op_input);
void OP_AddEnergy(struct s_polymer *p, int a1, int a2, double mul);
double Chi2(double *x, double *xexp, double *sigma, int n);
double OP_function(double **e, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms);
void OP_SamplePotential(struct s_polymer *p, struct s_mc_parms *parms, int ntypes, double **emat, int irun, double **ematold);
double OP_functionDiff(int iw, double deltae, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms);


