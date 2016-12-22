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

double Gauss(long *idum, double mean, double sigma)
{
 static int iset=0;
 static double gset;
 double fac, r, v1, v2;

 if (iset==0)
   {
    do {
        v1 = 2.0 * drand48() - 1.0;
        v2 = 2.0 * drand48() - 1.0;
                                       /** v1 = 2.0 * drand48() - 1.0;      
                                           v2 = 2.0 * drand48() - 1.0; **/
              r = v1*v1 + v2*v2;
       }
    while(r >= 1.0 || r == 0.0);
    fac = sqrt(-2.0*log(r)/r);
    gset = v1*fac;
    iset = 1;
    return(v2*fac*sigma + mean);
   }
 else
   {
    iset = 0;
    return (gset*sigma + mean);
   }

}


int Gaussian_Angles(double *angles,int n)
{
	int i;

      	for(i=0;i<n;i++)
	  angles[i]=Gauss((long*)42,0.0,1.0);	 	

	return 0;

}

void MatA( double **A, double **G, int dim,double a,double b)
{
     int i, j;
      
      for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
            { 
		if(i==j)
                	A[i][i]= 0.5*a*(1.+ b*G[i][j]);
                else
                	A[i][j]= 0.5*a*b*G[i][j]; 
              }
     return;
		   
		    
}


void Cholesky_2( double **L, double **A, int dim)
{
      int i,j,k;
      double sum=0.0;
      
      for( i=0; i<dim; i++)
            for( j=i; j<dim; j++)
            {

	          sum=A[i][j];
                  for(k=0;k<i;k++)
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

	for(i=0;i<natom_fragment;i++)
        	CopyBack((((p+ip)->back)+iw+i), (((f+ip)->back)+i));
       
}





int B_Metropolis(double deltaE,double temp,double WN,double WD,struct s_tables *t)
{
	if (deltaE<=0)
        	return 1;
			
	if(frand()<(WN/WD)*FastExp(-deltaE/temp,t))
	    	return 1;
      	

     return 0;
}




      



int CheckDistance(struct s_polymer *p,int ip,int iw,int nmul){
  
  double distance=Dist2( (((p+ip)->back)+(iw+nmul-1))->pos,(((p+ip)->back)+(iw+nmul))->pos );
  double d2next=( ( (p+ip)->back)+iw+nmul-1)->d2_next;

  if( DAbs(distance-d2next)<0.05)
    return 1;

  else return 0; 
  
}

int CheckAngles(struct s_polymer *p,int ip,int iw,int nmul){

  double angle;
  int out;


  angle=Angle( (((p+ip)->back)+iw+nmul-2)->pos, (((p+ip)->back)+iw+nmul-1)->pos, (((p+ip)->back)+iw+nmul)->pos,p->tables,&out);
 
  if(DAbs(angle-(((p+ip)->back)+iw+nmul-1)->a_next) > 2.0) return 0;
  
  
  angle=Angle( (((p+ip)->back)+iw+nmul-1)->pos, (((p+ip)->back)+iw+nmul)->pos, (((p+ip)->back)+iw+nmul+1)->pos,p->tables,&out);
  
  if ( DAbs(angle-(((p+ip)->back)+iw+nmul)->a_next) > 2.0) return 0;  

  else return 1; 
}

int CheckDihedral(struct s_polymer *p,int ip,int iw,int nmul){

  double dihedral;
  double dihcheck;
  int out;
  
  if((((p+ip)->back)+iw+nmul-1)->move==0){

    dihedral=Dihedral( (((p+ip)->back)+iw+nmul-3)->pos, (((p+ip)->back)+iw+nmul-2)->pos, (((p+ip)->back)+iw+nmul-1)->pos, (((p+ip)->back)+iw+nmul)->pos, p->tables, &out );
    dihcheck=(((p+ip)->back)+iw+nmul-1)->d_next;

  }
  else if((((p+ip)->back)+iw+nmul-2)->move==0){


    dihedral=Dihedral( (((p+ip)->back)+iw+nmul-1)->pos, (((p+ip)->back)+iw+nmul)->pos, (((p+ip)->back)+iw+nmul+1)->pos, (((p+ip)->back)+iw+nmul+2)->pos, p->tables, &out );
    dihcheck=(((p+ip)->back)+iw+nmul+1)->d_next;

  
    
  }
  else if((((p+ip)->back)+iw+nmul-3)->move==0){


     dihedral=Dihedral( (((p+ip)->back)+iw+nmul-2)->pos, (((p+ip)->back)+iw+nmul-1)->pos, (((p+ip)->back)+iw+nmul)->pos, (((p+ip)->back)+iw+nmul+1)->pos, p->tables, &out );
     dihcheck=(((p+ip)->back)+iw+nmul)->d_next;
    
  }
  else return 1;
  
  if( DAbs( DAbs(dihedral)-DAbs(dihcheck) ) > 2.0 ) return 0;
  else return 1;
  
}


/********************************************************************
 
 As PivotForward (geometry.c) :
 Pivot n atoms after iw of dih dw (moves from iw+2 on)
 affecting dihedral of iw+1
 
 pivots using initial coordinates of polymer p
 new coordinates on the fragment 
 p coordinates unchanged
 ********************************************************************/
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
	
	

        ifrag=0;
	for (i=iw+2;i<iw+2+n;++i)
	{

		U.x = (((p->back)+i)->pos).x - C.x;
		U.y = (((p->back)+i)->pos).y - C.y;
		U.z = (((p->back)+i)->pos).z - C.z;

		// rotate
		W.x = M[0][0] * U.x + M[1][0] * U.y + M[2][0] * U.z;
		W.y = M[0][1] * U.x + M[1][1] * U.y + M[2][1] * U.z;
		W.z = M[0][2] * U.x + M[1][2] * U.y + M[2][2] * U.z;

		
		(((fragment)->back)+ifrag)->pos.x = W.x + C.x;
		(((fragment)->back)+ifrag)->pos.y = W.y + C.y;
		(((fragment)->back)+ifrag)->pos.z = W.z + C.z;
   
		ifrag++;
	  
	}
	

	return 1;
}


int G_Matrix(struct s_polymer *fragment,struct s_polymer *p,int ip,int iw,int nmul,int nang,struct s_mc_parms *parms)
{
  int i,j,l;
  
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  double deriv1[nang][3],deriv2[nang][3],deriv3[nang][3];
  double fact,dw;
  
  
  x1=(((fragment+ip)->back)+nmul-3)->pos.x;
  y1=(((fragment+ip)->back)+nmul-3)->pos.y;
  z1=(((fragment+ip)->back)+nmul-3)->pos.z;
  
  x2=(((fragment+ip)->back)+nmul-2)->pos.x;
  y2=(((fragment+ip)->back)+nmul-2)->pos.y;
  z2=(((fragment+ip)->back)+nmul-2)->pos.z;
  
  x3=(((fragment+ip)->back)+nmul-1)->pos.x;
  y3=(((fragment+ip)->back)+nmul-1)->pos.y;
  z3=(((fragment+ip)->back)+nmul-1)->pos.z;
  
 
  

  dw=2.0;
  fact=1.0/dw;

  j=0;
  int k=0;
  for(i=iw;i<iw+nmul;++i)
  {
    
    
    
    k++;
    if ( (((p+ip)->back)+i)->move == 0 )  
      continue; 
     
    int nmoved=(nmul-2)-(k)+2;

    PivotForwardOnFragment((p+ip),(fragment+ip),i-1,dw,nmoved,parms);



    if ( 2 < nmoved)
    {
      deriv1[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-3-k)->pos.x )-x1)*fact;
      deriv1[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-3-k)->pos.y )-y1)*fact;
      deriv1[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-3-k)->pos.z )-z1)*fact;
    }
    
    else
    {
      deriv1[j][0]=0;
      deriv1[j][1]=0;
      deriv1[j][2]=0;           
    }
    
    if( 1 < nmoved)
    {
    deriv2[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+1)->pos.x )-x2)*fact;
    deriv2[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+1)->pos.y )-y2)*fact;
    deriv2[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+1)->pos.z )-z2)*fact;
    }

    else
    {
     deriv2[j][0]=0;
     deriv2[j][1]=0;
     deriv2[j][2]=0;
    }
    
    deriv3[j][0]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+2)->pos.x )-x3)*fact;
    deriv3[j][1]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+2)->pos.y )-y3)*fact;
    deriv3[j][2]=(((( (fragment+ip)->back)+parms->nmul_local-3-k+2)->pos.z )-z3)*fact;    
    
     
    j++;    
    
    
  }
  
   for(i=0;i<j;++i)
     for(k=0;k<j;++k)
       (fragment+ip)->G[i][k]=0.;
    
  
  for(i=0;i<j;++i)
    for(k=0;k<j;++k)
      for(l=0;l<3;++l)      			                  
        (fragment+ip)->G[i][k]+=(deriv1[i][l]*deriv1[k][l])+(deriv2[i][l]*deriv2[k][l])+(deriv3[i][l]*deriv3[k][l]);
  
  return 0; 
}





	  
int localpivot(int ip,int iw,struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
  
  
  double deltaE,W1=0.,W2=0.;

  int nang,i,m,ok=1;

  if (iw < nmul+3)
  {
    if ( (((p+ip)->back)+iw+1)->move == 0 ) { parms->mov++; return 0;}

    deltaE = -GetEnergyMonomerRange(p,0,iw,ip);							// two-body energy
  
    if (!parms->nodihpot)
      deltaE -= (((p+ip)->back)+iw+1)->e_dih;							// old dihedral energy
    if (!parms->nohfields)
      deltaE -= GetHFields(p,ip);										// old H fields energy
  
    ok = PivotBackward((p+ip),iw+1,0.2*(0.5-frand()),iw,parms);							// moved atoms are in [0,iw-1]; changed dihedral is iw+1

    ok *= AddSidechain(p,0,iw,ip);
    
    deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);			// new energy of iw with the others
    
    if (!parms->nodihpot)
      deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);
    if (!parms->nohfields)
	for(i=0;i<(p+ip)->nback;++i)
          deltaE += EnergyHFields(p,pot,parms,i,ip,1);

    ok *= Metropolis(deltaE,t,p->tables);

    if (ok == 0) 			// move rejected
    {
      UpdateMonomerRange(oldp,p,0,iw,ip,0);                         
      return 0;	
    }
    else
    {
      UpdateMonomerRange(p,oldp,0,iw,ip,0); 
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }
	
    return 0;
    
  }
  
  
  else if ( iw > ( (p+ip)->nback -3 -nmul   )  ) // do small forward
  { 
     if ( (((p+ip)->back)+iw)->move == 0 ) { parms->mov++; return 0; }

     deltaE=-GetEnergyMonomerRange(p,iw,(p+ip)->nback-1,ip);
     
     if (!parms->nodihpot)
        deltaE -= (((p+ip)->back)+iw)->e_dih;		
    if (!parms->nohfields)
       deltaE -= GetHFields(p,ip);

    ok = PivotForward((p+ip),iw-1,0.2*(0.5-frand()),(p+ip)->nback-iw-1,parms);			// moved atoms are in [iw+1,nback-1]; changed dihedral of iw
   
    ok *= AddSidechain(p,iw,(p+ip)->nback-1,ip);
    
    deltaE += EnergyMonomerRange(p,pot,iw,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
   
    if (!parms->nodihpot)
	deltaE += EnergyDihedrals(p,pot,iw,ip,1);
   
    if (!parms->nohfields)
    {
      for(i=0;i<(p+ip)->nback;++i)
      {
	deltaE += EnergyHFields(p,pot,parms,i,ip,1);
      }
      
    }

    ok *= Metropolis(deltaE,t,p->tables);
    
    if (ok == 0) 			// move rejected
    {
      UpdateMonomerRange(oldp,p,iw,(p+ip)->nback-1,ip,0);                           // ... if moved second half
      return 0;	
    }
    else
    {
      UpdateMonomerRange(p,oldp,iw,(p+ip)->nback-1,ip,0); 
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }

    
  } // end of small forward
  
  
  else {
    
    nang=0;
    
    for(i=0;i<nmul;i++)
      nang+=(((p+ip)->back)+(iw+i) )->move;

    CopyFragment(p,fragment,iw,parms->nmul_local+3,ip);
    
    
    deltaE=-GetEnergyMonomerRange(p,iw,iw+nmul,ip);
    
    if (!parms->nodihpot)
      for(i=0;i<nmul+3;i++)
        deltaE -= (((p+ip)->back)+iw+i)->e_dih;							// old dihedral energy
    
    if (!parms->nohfields)
      deltaE -= GetHFields(p,ip);
    
    Gaussian_Angles((fragment+ip)->g_ang,nang);        
    G_Matrix(fragment,p,ip,iw,nmul,nang,parms);    
    MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);    
    Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);    
    W1=exp(-Squared_n_Norma((fragment+ip)->g_ang,nang))*DetTriang((fragment+ip)->L,nang);
    InvertTriang((fragment+ip)->Y,(fragment+ip)->L,nang);
    TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang,nang);
        
    m=0;
    
    for(i=0;i<nmul-1;i++)
    {
      if( (((p+ip)->back)+iw+i+1)->move==0) continue;
      PivotForward((p+ip),iw+i,(fragment+ip)->d_ang[m],nmul-i-2,parms);
      m++;      
    }

    ok*=CheckDistance(p,ip,iw,nmul);
    ok*=CheckAngles(p,ip,iw,nmul);
    ok*=CheckDihedral(p,ip,iw,nmul);
    
    
    if(ok==1)
    {
      
      if(!parms->nosidechains)
       	ok *= AddSidechain(p,iw,iw+nmul,ip);
      
      deltaE += EnergyMonomerRange(p,pot,iw,iw+nmul,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
      
      if (!parms->nodihpot)
        for(i=0;i<nmul+3;i++)
      	  deltaE+=EnergyDihedrals(p,pot,iw+i,ip,1);
                      
      
      
      if (!parms->nohfields)
      	for(i=0;i<(p+ip)->nback;++i)
	  deltaE += EnergyHFields(p,pot,parms,i,ip,1);
	            
      
      G_Matrix(fragment,p,ip,iw,nmul,nang,parms);
      MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
      Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);    
      TransposedMatOnVect((fragment+ip)->L,(fragment+ip)->d_ang,(fragment+ip)->g_ang,nang);
      W2=exp(-Squared_n_Norma((fragment+ip)->g_ang,nang))*DetTriang((fragment+ip)->L,nang);
      ok*=B_Metropolis(deltaE,t,W2,W1,p->tables);
      
      if(ok==1)
      {
	UpdateMonomerRange(p,oldp,iw,iw+nmul,ip,parms->shell);  
	parms->acc++;
	p->etot += deltaE;
	return 1;
      }
      else
      {
	UpdateMonomerRange(oldp,p,iw,iw+nmul,ip,parms->shell);
	return 0;
      }
            
    }
                
    else
    {
      UpdateMonomerRange(oldp,p,iw,iw+nmul,ip,parms->shell);
      return 0;
    }
    
  }// end of local forward
  
  
  return 0; 
  
}
	  

int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
        int ip,iw,idir;
	int ok=1;
	
	
	ip   = irand(parms->npol);
	iw   = irand( (p+ip)->nback );
		
	ok=localpivot(ip,iw,p,oldp,fragment,pot,nmul,parms,t);  	

            	
	parms->mov++;
        
	return ok;

}	


int BackRub(struct s_polymer *p, struct s_polymer *oldp,struct s_potential *pot,struct s_mc_parms *parms,double t)
{

  int ok=1;
  int iw,ip;
  int i,out;

  double dw,deltaE;

  ip=irand( parms->npol);
  iw=3*irand((p+ip)->nback/3)+1;
  //iw=16;
  dw = 2.0 * (0.5 - frand());

  if(iw<3)
  {
    return ok;

  }
  if(iw>(p+ip)->nback-3)
  {
    return ok;
  }
  

  //fprintf(stdout,"backrub: iw = %d \t type : %s \n",iw,(((p+ip)->back)+iw)->type);
  deltaE = -GetEnergyMonomerRange(p,iw-3,iw+3,ip);                                                 // two-body energy

  if (!parms->nodihpot)
  {
      deltaE -= (((p+ip)->back)+iw-3)->e_dih;                                                   // old dihedral energy
      deltaE -= (((p+ip)->back)+iw-2)->e_dih;
      
      deltaE -= (((p+ip)->back)+iw+3)->e_dih;
      deltaE -= (((p+ip)->back)+iw+4)->e_dih;
  }
    if (!parms->nohfields)
      deltaE -= GetHFields(p,ip);                   

  ok=BackFlip(p,iw,dw);

  if( DAbs( Angle( (((p+ip)->back)+iw-4)->pos, (((p+ip)->back)+iw-3)->pos, (((p+ip)->back)+iw-2)->pos,p->tables,&out)-(((p+ip)->back)+iw-3)->a_next) > 5.0)
    ok=0;
  if(ok!=0) 
    if( DAbs( Angle( (((p+ip)->back)+iw+2)->pos, (((p+ip)->back)+iw+3)->pos, (((p+ip)->back)+iw+4)->pos,p->tables,&out)-(((p+ip)->back)+iw-3)->a_next) > 5.0)
      ok=0;


  if(ok==1)
  {

    if(!parms->nosidechains)
          ok *= AddSidechain(p,iw-3,iw+3,ip);
  
    deltaE += EnergyMonomerRange(p,pot,iw-3,iw+3,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
  
    if (!parms->nodihpot)
    {
      deltaE += EnergyDihedrals(p,pot,iw-3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw-2,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+4,ip,1);
    }
  
     if (!parms->nohfields)
       for(i=0;i<(p+ip)->nback;++i)
         deltaE += EnergyHFields(p,pot,parms,i,ip,1);
  
      ok *= Metropolis(deltaE,t,p->tables);


    }


    if (ok == 0)                        // move rejected
    {
      UpdateMonomerRange(oldp,p,iw-3,iw+3,ip,0);
      return 0;
    }

    else
    {
      UpdateMonomerRange(p,oldp,iw-3,iw+3,ip,0);
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }

  return ok;

}



int BackSideRub(struct s_polymer *p, struct s_polymer *oldp,struct s_potential *pot,struct s_mc_parms *parms,double t)
{

  int ok=1;
  int iw,ir,ip;
  int i,out;

  double dw,deltaE;

  ip=irand( parms->npol);
  iw=3*irand((p+ip)->nback/3)+1;
  //iw=16;
  dw = 2.0 * (0.5 - frand());

  if(iw<3)
  {
    return ok;

  }
  if(iw>(p+ip)->nback-3)
  {
    return ok;
  }


  //fprintf(stdout,"backrub: iw = %d \t type : %s \n",iw,(((p+ip)->back)+iw)->type);
  deltaE = -GetEnergyMonomerRange(p,iw-3,iw+3,ip);                                                 // two-body energy

  if (!parms->nodihpot)
  {
      deltaE -= (((p+ip)->back)+iw-3)->e_dih;                                                   // old dihedral energy
      deltaE -= (((p+ip)->back)+iw-2)->e_dih;

      deltaE -= (((p+ip)->back)+iw+3)->e_dih;
      deltaE -= (((p+ip)->back)+iw+4)->e_dih;
  }
    if (!parms->nohfields)
      deltaE -= GetHFields(p,ip);

  ok=BackFlip(p,iw,dw);

 if( DAbs( Angle( (((p+ip)->back)+iw-4)->pos, (((p+ip)->back)+iw-3)->pos, (((p+ip)->back)+iw-2)->pos,p->tables,&out)-(((p+ip)->back)+iw-3)->a_next) > 5.0)
    ok=0;
  if(ok!=0)
    if( DAbs( Angle( (((p+ip)->back)+iw+2)->pos, (((p+ip)->back)+iw+3)->pos, (((p+ip)->back)+iw+4)->pos,p->tables,&out)-(((p+ip)->back)+iw-3)->a_next) > 5.0)
      ok=0;


  if(ok==1)
  {

    ir=irand( (((p+ip)->back)+iw)->nrot );
    (((p+ip)->back)+iw)->irot=ir;
    //    ok = AddSidechain(p,iw,iw,ip);



    if(!parms->nosidechains)
          ok *= AddSidechain(p,iw-3,iw+3,ip);

    deltaE += EnergyMonomerRange(p,pot,iw-3,iw+3,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);

    if (!parms->nodihpot)
    {
      deltaE += EnergyDihedrals(p,pot,iw-3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw-2,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+4,ip,1);
    }

     if (!parms->nohfields)
       for(i=0;i<(p+ip)->nback;++i)
         deltaE += EnergyHFields(p,pot,parms,i,ip,1);

      ok *= Metropolis(deltaE,t,p->tables);


    }


    if (ok == 0)                        // move rejected
    {
      UpdateMonomerRange(oldp,p,iw-3,iw+3,ip,0);
      return 0;
    }

    else
    {
      UpdateMonomerRange(p,oldp,iw-3,iw+3,ip,0);
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }

  return ok;

}


