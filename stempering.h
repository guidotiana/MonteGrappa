/*
 * adjust_st.h
 *
 *  Created on: Nov 24, 2009
 *      Author: guido
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_multimin.h>

#define KB 		1.0

#include "do_mhistogram.h"	//must be after #include "define.h" (if it is included)

#define NTBIN		200		// # of bins to invert <E(T)> in AddNewTemperature
#define CONVERGENCE	1E-9	// convergence criterion for probability maximization
#define CONVER_WEAK 1E-7	// after NITERMAX, accept if this weaker condition is met
#define NITERMAX	5000	// maximum number of iterations       "
#define NHISTOMAX	150		// maximum number of histograms to remember
#define NEBINMAX	10001	// maximum number of energy bins
#define NRESTARTS	10		// number of last snapshots to restart

//#define NEEDMEMORY

struct st_restart
{
	double step;
	int ntemp;
	double *temp;
	double *g;
	double *lg;
};

struct st_stuff
{
	char method[15];
	int nm;					// shortening for method (1=st, 2=adaptive)
	int ntempmax;
	int ntemp;
	int nstep;			// attempts to change temperature every nstep
	int nstep_adj;			// adjust temperatures and weights every nstep_adj
	int nprint;
	int npre;				// equilibration state before everything
	double emin;
	double emax;
	double ebin;
	int nbin;
	char nftout[50];
	double k;				// new temperature put at k-times the standard deviation of the lowest
	double pthresh;			// probability threshold used to decrease the number of temperatures
	int debug;
	int keepall;			// keep all history of histograms
	char nfthe[50];			// output file for thermodynamics
	char nfdos[50];			// output file for density of states
	char nfdumb[50];		// output file for reliable dumb averages
	char nfdumb2[50];		// output file for current dumb averages
	char nfhisto[50];		// histogram file
	int nprintt;			// print thermodynamics every nprintt steps
	int paranoid;			// 0=give only warnings when something goes wrong, 1=stop
	double hthresh;			// threshold on fraction of time in a single temperature
	char anneal[15];		// how to set lower temperature (sigma/prob)
	double p_new;			// wished exchange probability
	int restart;			// if to restart
	double binthresh;		// threshold on (normalized) content of single bin
	double tonlyadd;		// below this temperature, only add temperatures
	int nkeepold;			// keep only last histograms
	int sum;				// 1=in keepall, sum together histograms with same temperature
	int removet;			// 0=do not remove temperatures, 99=any removal, ?=maximum decrease of ntemp in a shot
	double tstop;			// stop the simulation when this temperature is reached
	double tnorm;			// when this temperature is reached do nrmal stempering
	int nfail1;				// number of failures to fill holes with data from lg
	int nfail2;				// number of failures to fill holes with extrapolated free energies
	double phthresh;		// threshold on prob_up and prob_down to accept sampling
	int ignoreb;			// in mhistogram ignore temperatures which prevent equation solving
	double deltat;			// in mhistogram discard histograms which are closer than deltat in temperature
	int gmethod;			// how to calculate optimal weights g	

	double *temp;			// temperatures
	double *g;			// weigths
	int itemp;			// actual temperature
	double *prob_up;
	double *prob_down;
	int *counts;
	FILE *ftout;
	double **h;			// the histogram h[temp][e]
	double **htmp;			// temporary histogram h[temp][e] for normalization purposes
	double *f;			// log Z, only needed to mhistogram
	double *current_lg;		// current log of density of states
	double *reliable_lg;		// last reliable log of density of states
	double *reliable_t;		// 	"	temperatures
	double *reliable_g;		//	"	weights
	int nreliable_t;
	int count_st;			// counter to attempt temperature change
	int count_adj;			// counter to adjust temperatures/weights
	int count_print;		// counter to print temperatures
	int count_printt;		// counter to print thermodynamics
	int iter;			// current iteration of adjust_st
	int *oldt_iter;			// label each old histo with when it was collected
	int failure;			// how many consecutive failed runs
    	double energy;          	// potential energy of the system
    	int *binok;			// if to use a given energy bin (set in MHistogram)

	double **oldh;			// past histograms
	double *oldt;			// related temperatures
	int noldt;				// how many old histograms
	int noadjust;                   // if 1, do not adjust temperatures
    	double **out;			// thermodynamics output
    	double **boltzp;        	// boltzmann distribution
    	FILE *fpacc;
	double ttarget;			// below this temperature, set to this and make last run
	int ttarget_harvest;		// if 1, print trajectory
	int printpdb;
	FILE *pdbf;
	char pdbnf[500];
    	struct st_restart *st_restart;
};



// stempering.c
struct st_stuff *InitSTempering(void);
int STempering(double energy, double step, struct st_stuff *p);
int DoSTempering(double energy, struct st_stuff *p);
void OrderTemperatures(int ntemp, double *temp, double *g, double *prob_up,
		double *prob_down, double **h, int nbin, int debug, int *counts);
void PrintAverageEnergies(FILE *fout, double **h, double *temp, int ntemp, double emin, double ebin, int nbin);
void PrintHistogram(char *filename, double **h, double *t, int ntemp, int nbin, double emin, double ebin);
void Restart(struct st_stuff *p);
int FindClosestT(double t, int ntemp, double *temp);
void PrintStatistics(struct st_stuff *p, double step);

// adjust_st.c
void ManageRestart(struct st_stuff *p, double step);
int AdjustSTempering(struct st_stuff *p, double step);
int AddNewTemperature(int ntemp, double *temp, double *g, double *lg, double emin, double ebin,
		int nbin, int debug, double k, int paranoid, int *binok, int gmethod);
void OptimalWeights(int ntemp, double *t, double *g, double emin, double ebin, int nbin, double *lg, int debug, int *binok, int gmethod);
double EstimatedJumpProb (const gsl_vector *v, void *params);
double MaximizeTempProb(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double **boltzp, int *binok, int gmethod);
int ReduceTemperatures(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double pthresh, int ntmax, int removet, double **boltzp, int *binok, int gmethod);
int AddNewTemperature2(int ntemp, double *temp, double *g, double *lg, double emin, double ebin, int nbin,
		int debug, double p_new, int paranoid, double **boltzp, int *binok, int gmethod);
double EstimatedJumpSingleProbability(double t1, double t2, double *lg, double emin, double ebin, int nbin,
		int debug, double **boltzp, int *binok, int gmethod);
int CheckHistogramsAdjust(int *counts, double *temp, double *g, double *reliable_lg, int *ntemp, double emin, double ebin,
		int nbin, double hthresh, int debug, int ntmax, double **h, double *prob_up, double *prob_down, double pthresh, int *binok,
		int gmethod);
void NormalizeHistograms(double **h, double **htmp, int ntemp, int nbin, double binthresh, int gmethod);
int CheckHistograms(int *counts, double *temp, int ntemp, double htresh, int debug);
void AddHistogramPile(int ntemp, double *temp, double **h, double **oldh, double *oldt, int *oldt_iter, int *noldt, int nbin, int n, int iter, int debug, int sum);
int CheckHistogramsAdjust2(int *counts, double *temp, double *g, int *ntemp, double hthresh, int debug,
		double *reliable_t, double *reliable_g, int *nreliable_t, double *prob_up, double *prob_down, double phthresh);
double ExtrapolateLinear(double newt, double t1, double g1, double t2, double g2);
int CheckHistogramsAdjust3(int *counts, double *temp, double *g, int *ntemp, int nbin, double hthresh, int debug, int ntmax,
		double **h, double *prob_up, double *prob_down, double phthresh);
void FilterHistograms(double **h, double **hout, int ntemp, int nbin, double binthresh);
void CutHistograms(double **h, int ntemp, double emin, double ebin, int nbin, int *nbinnew, double *eminnew, int debug);
double CalculateG(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok, int method);

// do_mhistogram.c
double *MHistogram(double **h, double *t, int ntemp, int nbin, double emin, double ebin, int debug, double kb, double *,
		double *, int paranoid, int ignoreb, FILE *ff, double deltaT, int *binok);

// memory.c
double **AlloDoubleMat(int l, int m);
#ifdef NEEDMEMORY
double *AlloDouble(int l);
int *AlloInt(int l);
#endif
struct st_restart *AlloRestart(int ntemp, int nbin, int nres);
