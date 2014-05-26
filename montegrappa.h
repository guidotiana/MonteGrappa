/*
 * montegrappa.h
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <mpi.h>

#define NVER		0
#define NSUBVER	 	11	

#define NRESMAX		1000		// maximum number of backbone atoms
#define NATOMMAX	60000		// maximum number of atoms
#define NSIDEMAX	15		// maximum number of sidechain atoms for a given backbone atom
#define NROTMAX		40		// maximum number of rotamers allowed
#define NCONTMAX	200		// maximum number of other atoms in contact with a given one (in polymer structure)
#define NCONTMAX2	8000		// twice the number of other atoms in contact with a given one (in the moving structure, i.e. >>NCONTMAX)
#define NSHELLMAX	100		// maximum number of atoms in the shell of a given one
#define NMOVES		9		// number of types of allowed moves
#define NDIHFMAX	361		// maximum number of bins in tabled dihedral potential

#define NTEMPMAX	20
#define NAAMAX		1000
#define NCHAINMAX	25

// tables to calculate fast trigonometrics, exp, etc.
//#define TABLED_F			// activate it
					// square root (bin=0.0025)
#define FSQRT_L		80000		// # of elements in the table
#define FSQRT_BINRC	800.		// recirpocal of binning of the table
#define FSQRT_MAX	100		// = L / BIN
#define FSQRT_HBIN	0.000625	// half bin
#define FTRIG_L		72000		// # of elements in the table
#define FTRIG_BINRC	200.		// reciprocal of binning of the table
#define FTRIG_HBIN	0.0025		// half bin
#define FEXP_L		10000		// # of elements in the table
#define FEXP_BINRC	1000.		// reciprocal of binning of the table
#define FEXP_MAX	10		// = L / BIN
#define FEXP_HBIN	0.0005		// half bin
#define FACOS_L		40000		// # of elements in the table
#define FACOS_BINRC	20000.		// reciprocal of binning of the table
#define FACOS_HBIN	0.000025	// half bin

//LM_DELTAA per la mossa biased-gaussian
#define LM_DELTAA	0.0001
#define PARAM_A   10
#define PARAM_B   200
#define NMUL 8
#define NSTART 4
#define NCHECK 1
//#define DEBUG

#define PI 3.14159265
#define EPSILON 0.0001
#define EPSILON2 0.0001
#define LARGE 9E9

#define OPTIMIZEPOT


#include "struct_op.h"


#include "struct.h"

#include "optimizepot.h"


//#define STEMPERING


#define ACTIVE_MPI


#ifdef ACTIVE_MPI
#include "MPIfunc.h"
#define buffer_max 5000000 

#endif







void Welcome(FILE *fp);

// io.c
void Error(char *text);
void ReadPolymer(char *fname, struct s_polymer *p, FILE *flog, int npol, int debug, int *iamax, int *itypemax);
int FindKeyword(char *string, char *keyword);
void PrintPolymer(char *fname, struct s_polymer *p, int nchains);
void PrintPDBStream(struct s_polymer *p, int npol, FILE *fp);
struct s_mc_parms *ReadMcParms(char *fname);
int ReadPotential(char *fname, struct s_potential *u, struct s_mc_parms *parms, int na, int ntype);
void SetLookbackTables(struct s_polymer *p, int nc);
void PrintVector(FILE *fp, char c[], struct vector x);
void PrintAngles(FILE *fp, char c[], struct angles x);
void PrintStructure(struct s_polymer *p, int npol, FILE *fp, int shell);
void ReadParD(char *s, char key[20], int *par);
void ReadParL(char *s, char key[20], long *par);
void ReadParF(char *s, char key[20], double *par);
void ReadParS(char *s, char key[20], char *par);
void ReadParN(char *s, char key[20], int *par);
void PrintPotential(struct s_potential *u, char *eoutfile, int nat, int ntypes, int noangpot, int nodihpot, int hb);


// memory.c
struct s_polymer *AlloPolymer(int npol, int n, int nside, int nrot, int natoms, int shell, int noside, FILE *flog);
struct s_tables *InitTables(FILE *fp);
struct s_potential *AlloPotential(int natoms, int ntypes, int noangpot, int nodihpot, int hb);
int **AlloIntMatrix(int l, int m);
int *AlloInt(int l);
double **AlloDoubleMatrix(int l, int m);
double *AlloDouble(int l);

// mc.c
void Do_MC(struct s_polymer *p, struct s_polymer *fragment, struct s_polymer *replica, struct s_polymer *native, struct s_potential *pot, struct s_mc_parms *parms, FILE *ftrj, FILE *fe, struct s_polymer *oldp, FILE *fproc, int my_rank, int irun, MPI_Datatype Backtype, MPI_Datatype Sidetype, MPI_Datatype Rottype, MPI_Datatype Vecttype,MPI_Status astatus);
void CopyResiduePositions(struct s_back *from, struct s_back *to);
int MoveBackboneFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *u,int step, int debug, double t);
int MoveBackbonePivot(struct s_polymer *p,struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *mc_parms, double t);
int MoveMultiplePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *mc_parms, double t);
int Metropolis(double deltaE, double T, struct s_tables *t);
void UpdateMonomer(struct s_polymer *from, struct s_polymer *to, int w, int n, int shell);
int MoveSidechain(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *pot,
		int istep, int debug, double t);
void UpdateMonomerRange(struct s_polymer *from, struct s_polymer *to, int wfrom, int wto, int p, int shell);
int MoveLoosePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t);
void SoftExit();
int MoveMultipleFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		int istep, int debug, double t);
int MoveCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		int istep, int debug, double t);
void Anneal(struct s_mc_parms *p, double *t, int *counter, int *status, int *ok, int *ishell, int *mcount);

//bias_mp.c
void CopyResidueCoordinates(struct s_back *from,struct s_back *to);
int ComputeG(double **g,struct s_polymer *fragment,struct s_polymer *p,int k,int n,int npol,struct s_mc_parms *parms);
int MoveBiasedGaussian(struct s_polymer *p,struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *mc_parms,double t);
int MoveReallyLocal(struct s_polymer *p,struct s_polymer *oldp,struct s_potential *pot,int nmul,struct s_mc_parms *mc_parms,double t);
int CopyFragment(struct s_polymer *p,struct s_polymer *f,int k,int n,int ip);
int B_Metropolis(double deltaE,double T,double WN,double WD,struct s_tables *t);
struct s_polymer *AlloFragment(int npol,int nmul_local,FILE *flog);

//misc.c
int irand(int r);
double frand(void);
long Randomize(int n);
double FastSqrt(double x, struct s_tables *t);
double FastSin(double x, struct s_tables *t);
double FastCos(double x, struct s_tables *t);
double FastExp(double x, struct s_tables *t);
double FastAcos(double x, struct s_tables *t);
double Norm2(struct vector a);
void MatrixVectorProduct(double T[][3], struct vector in, struct vector *out);
double DAbs(double x);
int Abs(int x);
void CopyDoubleMatrix(double **from, double **to, int n, int m);
double gauss(double av, double sig);
void CopyVector(struct vector *from, struct vector *to);

// geometry.c
int Flip(struct s_polymer *p, int iw, double dw);
int MoveHead(struct s_polymer *p, struct s_mc_parms *mc_parms);
int MoveTail(struct s_polymer *p, struct s_mc_parms *mc_parms);
struct angles Cartesian2Spherical( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out);
struct vector Spherical2Cartesian(struct vector A, struct vector B,
                                  struct vector C, struct angles W, struct s_tables *tables, int *out);
struct vector RotateVector(struct vector v, double theta, int w, struct s_tables *tables);
int PivotForward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms);
int PivotBackward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms);
int Pivot(struct s_polymer *p, int iw, double dw, int ip, struct s_mc_parms *parms);
double Dist2(struct vector a, struct vector b);
double Dist(struct vector a, struct vector b);
double Angle( struct vector B, struct vector C, struct vector D, struct s_tables *tables, int *out);
double Dihedral( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out);
int AddSidechain(struct s_polymer *p, int istart, int istop, int ipol);
int FlipFragment(struct s_polymer *p, int ifrom, int ito, double dw);
void CopyPolymer(struct s_polymer *from, struct s_polymer *to, int cfrom, int cto, int noside, int norot);
void CopyAllPolymers(struct s_polymer *from, struct s_polymer *to, int n, int noside, int norot);
void DisplaceCoM(struct s_polymer *p, int ip, double dx, double dy, double dz);
double CosAngle( struct vector B, struct vector C, struct vector D, struct s_tables *tables, int *out);

// potential.c
double TotalEnergy(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int update, int sidechains, int debug, int iproc);
double EnergyMonomer(struct s_polymer *p, struct s_potential *u, int i, int ci, int npol, int update, int shell, int sidechains, int disentangle, int hb);
void AddContact(struct s_polymer *p, int i, int j, int ci, int cj, double e);
void AddShell(struct s_polymer *p, int i, int j, int ci, int cj);
void ResetContactsMonomer(struct s_polymer *p, int i, int ci);
double GetEnergyMonomer(struct s_polymer *p, int ip, int iw);
void PrintContacts(FILE *fp, struct s_polymer *p, int step);
double EnergyMonomerShell(struct s_polymer *p, struct s_potential *u, int i,int ip, int update);
void UpdateShell(struct s_polymer *p, struct s_mc_parms *parms);
double GetEnergyMonomerRange(struct s_polymer *p, int from, int to, int ip);
double EnergyMonomerRange(struct s_polymer *p, struct s_potential *u, int from, int to, int ip, int npol, int shell, int update, int nosidechains, int disentangle, int hb);
double EnergyPair(struct s_polymer *p, struct s_potential *u, int i, int j, int ci, int cj, int update, int sidechains, int disentangle, int tooclose, int hb);
double EnergyAngles(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update);
double EnergyDihedrals(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update);
void PrintEnergies(FILE *fp, int nc, int step, struct s_polymer *p);
void PrintEnergies_Parallel(FILE *fp, int nc, int step, struct s_polymer *p,int my_rank);
void CopyPotential(struct s_potential *from, struct s_potential *to, int nat, int ntypes);
int CheckOverlaps(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int nosidechains, int pr,FILE *fproc);
void CompareStructures(struct s_polymer *a, struct s_polymer *b, int nc, int natoms);
int EnergyBox(struct s_polymer *p, struct s_potential *u, int iw, int ic);
int EnergyBoxPolymer(struct s_polymer *p, struct s_potential *u, int ic);


// constants for gaussian random generator
float ran1(long *idum);
double gasdev(long *idum);

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
