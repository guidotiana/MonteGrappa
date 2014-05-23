/*
 * memory.c
 *
 *  Created on: Nov 18, 2009
 *      Author: guido
 */

#include "stempering.h"

/*****************************************************************************
 Allocate double matrix in the form x[l][m]
 *****************************************************************************/
double **AlloDoubleMat(int l, int m)
{
  double **x;
  int i,j;

  x = (double **) calloc(l,sizeof(double *));
  if (!x) { fprintf(stderr,"\n\nmatrix %dx%d",l,m); FatalError("Cannot allocate double matrix (1)");}

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(double));
      if (!(*(x+i))) { fprintf(stderr,"\n\nmatrix %dx%d",l,m); FatalError("Cannot allocate double matrix (2)");}
      for (j=0;j<m;j++) *(*(x+i)+j) = 0.;
    }

  return x;
}

#ifdef NEEDMEMORY
/*****************************************************************************
 Allocate double vector
 *****************************************************************************/
double *AlloDouble(int l)
{
  double *x;
  int i;

  x = (double *) calloc(l,sizeof(double));
  if (!x) { fprintf(stderr,"\n\nvector %d",l); FatalError("Cannot allocate double vector");}

  for (i=0;i<l;i++) *(x+i) =0;

  return x;
}

/*****************************************************************************
 Allocate int vector
 *****************************************************************************/
int *AlloInt(int l)
{
  int *x;
  int i;

  x = (int *) calloc(l,sizeof(int));
  if (!x) FatalError("Cannot allocate int");

  for (i=0;i<l;i++) *(x+i) =0;

  return x;
}

#endif

/*****************************************************************************
 Allocate structure for restart
 *****************************************************************************/
struct st_restart *AlloRestart(int ntemp, int nbin, int nres)
{
	struct st_restart *x;
	int i;

	x = (struct st_restart *) calloc(nres,sizeof(struct st_restart));
	if (!x) FatalError("Cannot allocate st_restart");

	for (i=0;i<nres;i++)
	{
		(x+i)->temp = calloc(ntemp,sizeof(double));
		if (!(x+i)->temp) FatalError("Cannot allocate st_restart (1)");

		(x+i)->g = calloc(ntemp,sizeof(double));
		if (!(x+i)->g) FatalError("Cannot allocate st_restart (2)");

		(x+i)->lg = calloc(nbin,sizeof(double));
		if (!(x+i)->lg) FatalError("Cannot allocate st_restart (3)");

		(x+i)->step = -1;
	}

	return x;
}

/*
int irand(int r)
{
  int i;
  double q;


  q=(double) rand()/(RAND_MAX+1.0);
  q *= r;
  i= (int) q;
  return i;
}

double frand()
{
  double q;

  q=(double) rand()/(RAND_MAX+1.0);
  return q;
}
*/
