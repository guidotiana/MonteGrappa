/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
    Copyright (C) 2014  G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, R. Meloni

    MonteGrappa is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    This program is free. We just ask you to cite in the publications you get 
    making use of    MonteGrappa the following article:

    G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, E. Heard,
    L. Giorgetti and R. Meloni, "MonteGrappa: an iterative Monte Carlo program
    to optimize biomolecular potentials in simplified models" Comp. Phys. Comm.
    Volume 186, January 2015, Pages 93–104.     DOI: 10.1016/j.cpc.2014.09.012

   
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "montegrappa.h"



struct s_polymer *Allo_Fragment(int npol, int nback,int nangmax, FILE *flog)
{
	
	int ipol;
	struct s_polymer *f;

	f=(struct s_polymer *)calloc(npol,sizeof(struct s_polymer));
	if(!f)      Error("Cannot Allocate fragment");


	for(ipol=0;ipol<npol;++ipol)
	{
		
		(f+ipol)->nback=nback;
		(f+ipol)->back=(struct s_back *)calloc((nback+3),sizeof(struct s_back));
		if(!(f+ipol)->back)     Error("\tAlloFragment(): Cannot Allocate backbone");

                (f+ipol)->A=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->G=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->L=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->Y=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->g_ang=AlloDouble(nangmax);      
                (f+ipol)->d_ang=AlloDouble(nangmax);
	
	}

	return f;

}

double DetTriang( double **Mat, int n )
{

	int i;
	double det=1.0;
	
	for(i=0; i<n; ++i)
   		det = det * Mat[i][i];

	return det;

}


void FreeFragment(struct s_polymer *fragment,struct s_mc_parms *parms)
{

	int i;
	for(i=0;i<parms->npol;++i)
        {
	        FreeDoubleMatrix((fragment+i)->A,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->G,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->L,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->Y,parms->nmul_local);
                free((fragment+i)->d_ang );
                free((fragment+i)->g_ang);
		free((fragment+i)->back);	
	}

	free(fragment);

}


void InvertTriang( double **Inv, double **Mat, int n )
{
	int i,k,a;
	
	for(i=0; i<n; ++i)
   		for(k=0; k<n; ++k)
     			Inv[i][k]=0.0;

	for(i=0; i<n; ++i)
   		Inv[i][i] = 1./Mat[i][i];

	for(i=0; i<n; ++i)
  		for(a=1; a<=i; a++)
  		{
     			for( k= i-a; k<i; ++k)
       				Inv[i][i-a] += (Mat[i][k] * Inv[k][i-a]);
       			Inv[i][i-a] /= -Mat[i][i];
  		}

	return;

}




double Squared_n_Norma( double *vect, int dim )
{
	int i;
	double norm=0.0;

	for(i=0; i<dim; ++i)
    		norm += vect[i]*vect[i];
	
	return norm;
}

void TransposedMatOnVect( double **Mat, double *vect, double *risu, int n)
{
       int i,j;

      for(i=0; i<n; ++i)
      {
            risu[i] = 0.0;
            for(j=0; j<n; ++j)
                  risu[i] += Mat[j][i] * vect[j];
      }

      return;
}

int Gaussian_Angles(double *angles,int n)
{
	int i;
        double r1,r2;

      	for(i=0;i<n;++i)
        {
                r1=(double)rand()/(double)RAND_MAX;
                r2=(double)rand()/(double)RAND_MAX;
                angles[i]=sqrt(-log(r1))*cos(2*PI*r2);
      	}

     	return 0;

}

void MatA( double **A, double **G, int dim,double a,double b)
{

      int i, j;
      //#pragma omp parallel for private(j)
      for(i=0; i<dim; ++i)
            for(j=0; j<dim; ++j)
            {
                   
		if(i==j)
                	A[i][i]= a*(1+ b*G[i][j])/2;
                else
                	A[i][j]= a*b*G[i][j]/2;
            
              }
                    return;
}


void Cholesky_2( double **L, double **A, int dim)
{
      int i,j,k;
      double sum;
      
      for( i=0; i<dim; ++i)
            for( j=i; j<dim; ++j)
            {
                  sum=A[i][j];
                  for(k=0;k<i;++k)
                        sum += -L[i][k]*L[j][k];
                  if(i == j)
                  {
                        if(sum<=0.0)
                              Error("Cholesky decomposition failed (check the BGS parameters!)");
                        L[i][i]= sqrt(sum);
                  }
                  else
                        L[j][i]= sum/L[i][i];
            }

      return;
}
      
void CopyBack(struct s_back *from, struct s_back *to)
{
      (to)->nside=(from)->nside;
      
      (to->pos).x = (from->pos).x;
      (to->pos).y = (from->pos).y;
      (to->pos).z = (from->pos).z;

      (to)->move=(from)->move;

      strncpy(to->type,from->type,5);

  
}


void CopyFragment(struct s_polymer *p,struct s_polymer *f,int iw,int natom_fragment,int ip)
{
      
        int i;

	for(i=0;i<natom_fragment;++i)
        	CopyBack((((p+ip)->back)+iw+i), (((f+ip)->back)+i));
        
}





int B_Metropolis(double deltaE,double temp,double WN,double WD,struct s_tables *t)
{
	if (deltaE<=0)
        	return 1;
	
	double p;
	p=FastExp(-deltaE/temp,t);
	p=p*(WN/WD);
    
	double random=frand();

	if(random<p)
	{
        	return 1;
      	}

     return 0;
}




      
int Compute_G(struct s_polymer *fragment,struct s_polymer *p,int ip,int iw,int nmul,int nang,struct s_mc_parms *parms)
{

  /*
    int ok=1;
	int natom_fragment=nmul+3;
        int i,j,l,i_ang;


      	int k1;
      	int k2;
      	int k3;
      	double deriv1[nang][3],deriv2[nang][3],deriv3[nang][3];
      	double x1,y1,z1,x2,y2,z2,x3,y3,z3;
      	double fact=(1./(LM_DELTAA));


      	CopyFragment(p,fragment,iw,natom_fragment,ip);


      	k1=nmul-2;
      	k2=nmul-1;
      	k3=nmul;


      	x1=(((fragment+ip)->back)+(k1))->pos.x;
      	y1=(((fragment+ip)->back)+(k1))->pos.y;
	z1=(((fragment+ip)->back)+(k1))->pos.z;
      	x2=(((fragment+ip)->back)+(k2))->pos.x;
      	y2=(((fragment+ip)->back)+(k2))->pos.y;
      	z2=(((fragment+ip)->back)+(k2))->pos.z;
      	x3=(((fragment+ip)->back)+(k3))->pos.x;
      	y3=(((fragment+ip)->back)+(k3))->pos.y;
      	z3=(((fragment+ip)->back)+(k3))->pos.z;
      	i_ang=0;
      

      	for(i=0;i<natom_fragment;++i)
      	{
        	if( (((fragment+ip)->back)+i+1)->move==1 )
            	{
                  PivotForward((fragment+ip),i,LM_DELTAA,natom_fragment-i-2,parms);

                  deriv1[i_ang][0]=(( (((fragment+ip)->back)+(k1))->pos.x )-x1);
                  deriv1[i_ang][1]=( (((fragment+ip)->back)+(k1))->pos.y )-y1;
                  deriv1[i_ang][2]=( (((fragment+ip)->back)+(k1))->pos.z )-z1;
                  deriv2[i_ang][0]=( (((fragment+ip)->back)+(k2))->pos.x )-x2;
                  deriv2[i_ang][1]=( (((fragment+ip)->back)+(k2))->pos.y )-y2;
                  deriv2[i_ang][2]=( (((fragment+ip)->back)+(k2))->pos.z )-z2;
                  deriv3[i_ang][0]=( (((fragment+ip)->back)+(k3))->pos.x )-x3;
                  deriv3[i_ang][1]=( (((fragment+ip)->back)+(k3))->pos.y )-y3;
                  deriv3[i_ang][2]=( (((fragment+ip)->back)+(k3))->pos.z )-z3;

                  deriv1[i_ang][0]=(deriv1[i_ang][0]*fact);
                  deriv1[i_ang][1]=(deriv1[i_ang][1]*fact);
                  deriv1[i_ang][2]=(deriv1[i_ang][2]*fact);
                  deriv2[i_ang][0]=(deriv2[i_ang][0]*fact);
                  deriv2[i_ang][1]=(deriv2[i_ang][1]*fact);
                  deriv2[i_ang][2]=(deriv2[i_ang][2]*fact);
                  deriv3[i_ang][0]=(deriv3[i_ang][0]*fact);
                  deriv3[i_ang][1]=(deriv3[i_ang][1]*fact);
                  deriv3[i_ang][2]=(deriv3[i_ang][2]*fact);

                 CopyFragment(p,fragment,iw,nmul,natom_fragment,ip);

                  i_ang++;
                  if(i_ang==(nang)) break;
            }
      }



     for(i=0;i<nang;++i)
            for(j=0;j<nang;++j)
		(fragment+ip)->G[i][j]=0;

      for(i=0;i<nang;++i)
            for(j=0;j<nang;++j)
                  for(l=0;l<3;++l)
      			{                  
			(fragment+ip)->G[i][j]+=(deriv1[i][l]*deriv1[j][l])+(deriv2[i][l]*deriv2[j][l])+(deriv3[i][l]*deriv3[j][l]);

			}


//	fprintf(stderr,"\n\n\n\n");
//	for(i=0;i<nang;++i)
//	{
//	fprintf(stderr,"\n");
//	for(j=0;j<nang;++j)
//		fprintf(stderr,"%lf\t",(fragment+ip)->G[i][j]);

//	}

*/
      return 0;
      
}


int PivotForwardOnFragment(struct s_polymer *p,struct s_polymer *fragment, int iw, double dw, int n, struct s_mc_parms *parms)
{
	
	int i,ifrag;
	struct vector C,V,U,W;
	double norm,inv_norm,M[3][3],cw,sw,ocw;

	#ifdef DEBUG
	if (iw+n+1 >= p->nback) Error("attempting to pivot beyond the end of the chain");
	#endif

	// the pivoting atom
	C.x = (((p->back)+iw)->pos).x;
	C.y = (((p->back)+iw)->pos).y;
	C.z = (((p->back)+iw)->pos).z;

	// the normalized vector around which to pivot
	V.x = (((p->back)+iw+1)->pos).x - C.x;
	V.y = (((p->back)+iw+1)->pos).y - C.y;
	V.z = (((p->back)+iw+1)->pos).z - C.z;

	norm = FastSqrt(Norm2(V),p->tables);
	if (norm<EPSILON) return 0; 		//PARANOID
	inv_norm = 1. / norm;
	if (inv_norm<EPSILON) return 0; 	//PARANOID
	V.x *= inv_norm;
	V.y *= inv_norm;
	V.z *= inv_norm;


	// the rotation matrix
	cw = FastCos(dw,p->tables);
	sw = FastSin(dw,p->tables);
	ocw = 1. - cw;

	M[0][0] = V.x * V.x * ocw + cw;
	M[0][1] = V.x * V.y * ocw + V.z * sw;
	M[0][2] = V.x * V.z * ocw - V.y * sw;

	M[1][0] = V.x * V.y * ocw - V.z * sw;
	M[1][1] = V.y * V.y * ocw + cw;
	M[1][2] = V.y * V.z * ocw + V.x * sw;

	M[2][0] = V.x * V.z * ocw + V.y * sw;
	M[2][1] = V.y * V.z * ocw - V.x * sw;
	M[2][2] = V.z * V.z * ocw + cw;
	
	
	// move n atoms from iw+2 to iw+2+n

        ifrag=0;
	fprintf(stdout,"pivot forward on fragment, move %d atoms from %d to %d \n",n,iw+2,iw+1+n);
	for (i=iw+2;i<iw+2+n;++i)
     // for(i=iw;i<iw+n;++i)
	{
//		fprintf(stderr,"\ni = %d \n",i);
            // shift each to be centred in C
		U.x = (((p->back)+i)->pos).x - C.x;
		U.y = (((p->back)+i)->pos).y - C.y;
		U.z = (((p->back)+i)->pos).z - C.z;

		// rotate
		W.x = M[0][0] * U.x + M[1][0] * U.y + M[2][0] * U.z;
		W.y = M[0][1] * U.x + M[1][1] * U.y + M[2][1] * U.z;
		W.z = M[0][2] * U.x + M[1][2] * U.y + M[2][2] * U.z;

		// shift back
		fprintf(stdout,"DEBUG shift back backbone atom %d of fragment\n",ifrag);
		(((fragment)->back)+ifrag)->pos.x = W.x + C.x;
		(((fragment)->back)+ifrag)->pos.y = W.y + C.y;
		(((fragment)->back)+ifrag)->pos.z = W.z + C.z;
   
		ifrag++;
	  
	}
	for (i=iw+2;i<iw+2+n;++i)
	{
	fprintf(stdout,"DEBUG polymer coordinates backbone %d\t%lf\t%f\t%f\n",i,(((p)->back)+i)->pos.x,(((p)->back)+i)->pos.y,(((p)->back)+i)->pos.z);
	}
	
	
	for(i=0;i<ifrag;++i)
	{
	fprintf(stdout,"DEBUG fragment coordinates backbone %d\t%lf\t%f\t%f\n",i,(((fragment)->back)+i)->pos.x,(((fragment)->back)+i)->pos.y,(((fragment)->back)+i)->pos.z);
	}

	return 1;
}


int G_Matrix(struct s_polymer *fragment,struct s_polymer *p,int ip,int iw,int nmul,int nang,struct s_mc_parms *parms)
{
  int i,j,l;
  
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  double deriv1[nang][3],deriv2[nang][3],deriv3[nang][3];
  double fact;
  
  for(i=iw;i<iw+parms->nmul_local;++i)
  {
    fprintf(stdout,"DEBUG before copy polymer coordinates backbone %d\t%lf\t%f\t%f\n",i,(((p+ip)->back)+i)->pos.x,(((p+ip)->back)+i)->pos.y,(((p+ip)->back)+i)->pos.z);
  }
  
  CopyFragment(p,fragment,iw,parms->nmul_local,ip);
  
  
  for(i=0;i<parms->nmul_local;++i)
  {
    fprintf(stdout,"DEBUG after copy fragment coordinates backbone %d\t%lf\t%f\t%f\n",i,(((fragment+ip)->back)+i)->pos.x,(((fragment+ip)->back)+i)->pos.y,(((fragment+ip)->back)+i)->pos.z);
  }
  
  x1=(((fragment+ip)->back)+parms->nmul_local-3)->pos.x;
  y1=(((fragment+ip)->back)+parms->nmul_local-3)->pos.y;
  z1=(((fragment+ip)->back)+parms->nmul_local-3)->pos.z;
  
  x2=(((fragment+ip)->back)+parms->nmul_local-2)->pos.x;
  y2=(((fragment+ip)->back)+parms->nmul_local-2)->pos.y;
  z2=(((fragment+ip)->back)+parms->nmul_local-2)->pos.z;
  
  x3=(((fragment+ip)->back)+parms->nmul_local-1)->pos.x;
  y3=(((fragment+ip)->back)+parms->nmul_local-1)->pos.y;
  z3=(((fragment+ip)->back)+parms->nmul_local-1)->pos.z;
  
  fprintf(stdout,"DEBUG last x y z = %lf %lf %lf\n",x1,y1,z1);
  fprintf(stdout,"DEBUG last x y z = %lf %lf %lf\n",x2,y2,z2);
  fprintf(stdout,"DEBUG last x y z = %lf %lf %lf\n",x3,y3,z3);

  double dw;
  dw=parms->dw_pivot * (0.5 - frand());
  fact=1./dw;
  fprintf(stdout,"dw=%lf\n",dw);
  j=0;
  
  for(i=iw;i<iw+nmul;++i)
  {
    if ( (((p+ip)->back)+i)->move == 0 ) continue;
    PivotForwardOnFragment((p+ip),(fragment+ip),i,dw,nmul-2-j,parms);
    deriv1[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-3-j)->pos.x )-x1)*fact;
    deriv1[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-3-j)->pos.x )-y1)*fact;
    deriv1[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-3-j)->pos.x )-z1)*fact;
    deriv2[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-2-j)->pos.x )-x2)*fact;
    deriv2[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-2-j)->pos.x )-y1)*fact;
    deriv2[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-2-j)->pos.x )-z2)*fact;
    deriv3[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-1-j)->pos.x )-x3)*fact;
    deriv3[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-1-j)->pos.x )-y3)*fact;
    deriv3[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-1-j)->pos.x )-z3)*fact;    
        
    ++j;
    
  }
  
  for(i=0;i<j;++i)
  {
     fprintf(stdout,"DEBUG deriv1 %d = %lf %lf %lf\n",i,deriv1[i][0],deriv1[i][1],deriv1[i][2]);
     fprintf(stdout,"DEBUG deriv2 %d = %lf %lf %lf\n",i,deriv2[i][0],deriv2[i][1],deriv2[i][2]);
     fprintf(stdout,"DEBUG deriv3 %d = %lf %lf %lf\n",i,deriv3[i][0],deriv3[i][1],deriv3[i][2]);

    
     
  }
  fprintf(stdout,"DEBUG we performed %d pivots\n",j);

  for(i=0;i<nang;++i)
    for(j=0;j<nang;++j)
      for(l=0;l<3;++l)      			                  
        (fragment+ip)->G[i][j]+=(deriv1[i][l]*deriv1[j][l])+(deriv2[i][l]*deriv2[j][l])+(deriv3[i][l]*deriv3[j][l]);

			
  
  
 return 0; 
}




int localbackward(int ip,int iw,struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
  
  double deltaE;
  int nang,i;
  
  if ( iw < nmul)
  {
    
    // just pivot backward with small angle
    fprintf(stdout,"small backward\n");
    fprintf(stdout,"DEBUG get energy between 0 and %d\n",iw);
    deltaE=-GetEnergyMonomerRange(p,0,iw,ip);
    fprintf(stdout,"DEBUG energy is %lf\n",deltaE);

    nang=0;
    
    for(i=0;i<iw;++i)
    {
        fprintf(stdout,"get move of backbone %d=%d\n",iw-i-1,(((p+ip)->back)+(iw-i-1))->move);
	nang+=(((p+ip)->back)+(iw-i-1))->move;
    }

    fprintf(stdout,"DEBUG nang = %d\n",nang);
    
    
  }
  
  else
  {
   // do a backward local move 
        fprintf(stdout,"DEBUG local backward\n");
	fprintf(stdout,"DEBUG get energy between %d and %d\n",iw-nmul,iw);
    	deltaE=-GetEnergyMonomerRange(p,iw-nmul,iw,ip);
	fprintf(stdout,"DEBUG energy is %lf\n",deltaE);

	nang=0;
    
        for(i=0;i<nmul;++i)
	{
                fprintf(stdout,"get move of backbone %d=%d\n",iw-i-1,(((p+ip)->back)+(iw-i-1))->move);
		nang+=(((p+ip)->back)+(iw-i-1))->move;
	}

        fprintf(stdout,"DEBUG nang = %d\n",nang);  
    
  }
  
  return 0; 
}
	  
int localforward(int ip,int iw,struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
  double deltaE,e,psisquared,detA1,detL1,W1,rc2;
  rc2=parms->r_cloose*parms->r_cloose;
  int nang,i,m,j;
  
  if ( iw > ( (p+ip)->nback - nmul)-1 )
  {
     // just pivot forward with small angle
      fprintf(stdout,"small forward\n");
      fprintf(stdout,"DEBUG get energy between %d and %d\n",iw,(p+ip)->nback-1);
      deltaE=-GetEnergyMonomerRange(p,iw,(p+ip)->nback-1,ip);
      fprintf(stdout,"DEBUG energy is %lf\n",deltaE);

      nang=0;
   
      for(i=iw;i<((p+ip)->nback)-1;++i)
      {
                fprintf(stdout,"get move of backbone %d=%d\n",i+1,(((p+ip)->back)+(i+1))->move);

		nang+= (((p+ip)->back)+(i+1))->move;
	}

    fprintf(stdout,"DEBUG nang = %d\n",nang);  

    
  }
  
  else{
    // do a forward local move
    
    fprintf(stdout,"DEBUG local forward\n");
    fprintf(stdout,"DEBUG get energy between %d and %d\n",iw,iw+nmul);

    deltaE=-GetEnergyMonomerRange(p,iw,iw+nmul,ip);
    if(!parms->nodihpot)
      fprintf(stdout,"DEBUG QUA BISOGNA PRENDERE L'ENERGIA DIEDRI E ANGOLI DA METTERE\n");
    
    
    fprintf(stdout,"DEBUG energy is %lf\n",deltaE);
 
    nang=0;
    
    for(i=0;i<nmul;++i)
	{
                fprintf(stdout,"get move of backbone %d = %d\n",iw+i+1,(((p+ip)->back)+(iw+i+1) )->move);
      
		nang+=(((p+ip)->back)+(iw+i+1) )->move;
	}

    fprintf(stdout,"DEBUG nang = %d\n",nang);  
    
    Gaussian_Angles((fragment+ip)->g_ang,nang);
    fprintf(stdout,"\nDEBUG generated %d gaussian angles\n",nang);
    
    for(i=0;i<nang;++i)
    fprintf(stdout,"DEBUG   %d = %lf\n",i,(((fragment+ip)->g_ang[i])));
 
    fprintf(stdout,"DEBUG now compute G\n");
    
    G_Matrix(fragment,p,ip,iw,nmul,nang,parms);
  
    MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
    Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
    psisquared=Squared_n_Norma((fragment+ip)->g_ang,nang);
    e=exp(-psisquared);
    detL1=DetTriang((fragment+ip)->L,nang);
    detA1=detL1*detL1;
    W1=e * sqrt(detA1);
    InvertTriang((fragment+ip)->Y,(fragment+ip)->L,nang);
    TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang,nang);

    fprintf(stdout,"DEBUG pivoting polymer using d-angles\n");
    m=0;
    for(i=0;i<nmul;++i)
    {
      
      if( (((p+ip)->back)+iw+i+1)->move==0) continue;
      fprintf(stdout,"DEBUG pivoting after %d, using angle %d\n",iw+i,m);
      for(j=iw;j<iw+nmul;++j)
      {
        fprintf(stdout,"DEBUG   before pivoting polymer coordinates backbone %d\t%lf\t%f\t%f\n",j,(((p+ip)->back)+j)->pos.x,(((p+ip)->back)+j)->pos.y,(((p+ip)->back)+j)->pos.z);
      }
      PivotForward((p+ip),iw+i,(fragment+ip)->d_ang[m],nmul-i-2,parms);
      for(j=iw;j<iw+nmul;++j)
      {
        fprintf(stdout,"DEBUG   after pivoting polymer coordinates backbone %d\t%lf\t%f\t%f\n",j,(((p+ip)->back)+j)->pos.x,(((p+ip)->back)+j)->pos.y,(((p+ip)->back)+j)->pos.z);
      }
      ++m;
      //if(m==nang)
      //	break;
      
	
    }
    
    fprintf(stdout,"DEBUG end of pivots. %d angles were used\n",m);
   
    fprintf(stdout,"DEBUG now check distance between moved segment and fixed part\n");
    fprintf(stdout,"DEBUG dist between %d and %d\n",iw+nmul-1,iw+nmul);
    /*
    if(DAbs (Dist2( (((p+ip)->back)+iw+nmul)->pos,(((p+ip)->back)+iw+nmul+1)->pos )) -(((p+ip)->back)+iw)->d2_next  < rc2)
    {
     fprintf(stdout,"DEBUG distance is not good! %lf > %lf\n",(Dist2( (((p+ip)->back)+iw)->pos,(((p+ip)->back)+iw+1)->pos ) -(((p+ip)->back)+iw),rc2); 
    }
    */
    
  }
    
  

  
     
  
  return 0; 
}
	  

int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
        int ip,iw,idir;
	int ok=1;
	
	
	ip   = irand(parms->npol);
	iw   = irand( (p+ip)->nback );
        idir = rand() % 2 ? 1 : -1;
	
	
	fprintf(stdout,"DEBUG selected polymer %d\n",ip);
        fprintf(stdout,"DEBUG selected backbone atom %d\n",iw);
        fprintf(stdout,"DEBUG selected direction %d\n",idir);
	
	if(idir==-1)
	  ok=localbackward(ip,iw,p,oldp,fragment,pot,nmul,parms,t);
	else
	  ok=localforward(ip,iw,p,oldp,fragment,pot,nmul,parms,t);
	  	
/*	
	
	deltaE = -GetEnergyMonomerRange(p,iw-nmul+1,iw+1,ip);
  

		if((((p+ip)->back)+iw+1)->iapdb==0)
			iapdbtocheck=1;
		if((((p+ip)->back)+iw+1)->iapdb==1)
			iapdbtocheck=2;
		if((((p+ip)->back)+iw+1)->iapdb==2)
			iapdbtocheck=0;	
		int out;
		Gaussian_Angles((fragment+ip)->g_ang,nang);

		Compute_G(fragment,p,ip,iw-nmul+1,nmul,nang,parms);	
		MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
		Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
		psisquared=Squared_n_Norma((fragment+ip)->g_ang,nang);
                e=exp(-psisquared);
		detL1=DetTriang((fragment+ip)->L,nang);
		detA1=detL1*detL1;
		W1=e * sqrt(detA1);
		InvertTriang((fragment+ip)->Y,(fragment+ip)->L,nang);
		TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang,nang);

		if(!parms->nodihpot)
			for(i=0;i<nmul+3;++i)
			{
				deltaE -= (((p+ip)->back)+iw-nmul+i+1)->e_dih;
			}

			if(!parms->noangpot)
			{		
		
                	deltaE-=(((p+ip)->back)+iw)->e_ang;
  	            	deltaE-=(((p+ip)->back)+iw+1)->e_ang;
			}


		m=0;
		for(i=0;i<nmul;++i)
		{
			if( (((p+ip)->back)+iw+i-nmul+2)->move==1)
			{
				ok*=PivotForward((p+ip),iw-nmul+i+1,(fragment+ip)->d_ang[m],nmul-i-2,parms);
				++m;
				if(m==nang)
				break;
                  	}
	
		}

                                                                   

		
		double rc2=parms->r_cloose*parms->r_cloose;
		double dihedral=Dihedral( (((p+ip)->back)+iw-iapdbtocheck)->pos, (((p+ip)->back)+iw+1-iapdbtocheck)->pos, (((p+ip)->back)+iw+2-iapdbtocheck)->pos, (((p+ip)->back)+iw+3-iapdbtocheck)->pos, p->tables, &out );

		
		if( DAbs (Dist2( (((p+ip)->back)+iw)->pos,(((p+ip)->back)+iw+1)->pos ) -(((p+ip)->back)+iw)->d2_next < rc2 )  &&   DAbs( Angle( (((p+ip)->back)+iw-1)->pos, (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (p+ip)->tables, &out)  - (((p+ip)->back)+iw-1)->a_next) < 0.2 && DAbs(180-DAbs(dihedral)) < DELTAOMEGA)
		{		

			if(!parms->nosidechains)
			{
				ok *= AddSidechain(p,iw-nmul+1,iw+1,ip);
			}
	
			deltaE += EnergyMonomerRange(p,pot,iw-nmul+1,iw+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);

			if (!parms->nodihpot)
        	        	for(i=0;i<nmul+3;++i)
				{
	                      		deltaE+=EnergyDihedrals(p,pot,iw-nmul+i+1,ip,1);
				}                  

                  	if(!parms->noangpot)
                  	{
				//		fprintf(stderr,"\nCOMPUTING ANGLE ENERGY deltaE=%lf\t",deltaE);
        			deltaE+=EnergyAngles(p,pot,iw,ip,1);
                       		deltaE+=EnergyAngles(p,pot,iw+1,ip,1);
				//		fprintf(stderr,"-> %lf\n",deltaE);
 	              	}

       		
//			fprintf(stderr,"dietro\n");
                  	Compute_G(fragment,p,ip,iw-nmul+1,nmul,nang,parms);
 	          	MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
			Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
			detL2=DetTriang((fragment+ip)->L,nang);
			detA2=detL2*detL2;
			TransposedMatOnVect((fragment+ip)->L,(fragment+ip)->d_ang,(fragment+ip)->g_ang,nang);
			psisquared=Squared_n_Norma((fragment+ip)->g_ang,nang);
			e=exp(-psisquared);
			W2=e*sqrt(detA2);
			ok*=B_Metropolis(deltaE,t,W2,W1,p->tables);

			if(ok==0) //move rejected
                  	{
      				UpdateMonomerRange(oldp,p,iw-nmul+1,iw+1,ip,parms->shell);
			}
                  	else	//move accepted
                  	{

                	        UpdateMonomerRange(p,oldp,iw-nmul+1,iw+1,ip,parms->shell);
                		
			        p->etot+=deltaE;
       			
		                parms->acc++;
                  	}
			
			
		}//end of loose condition

		else //if not loose pivot
        	{
			ok=0;
	        	UpdateMonomerRange(oldp,p,iw-nmul+1,iw+1,ip,parms->shell);
            	}
*/
            	
	parms->mov++;

	return ok;

}	

