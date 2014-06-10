#include "montegrappa.h"



struct s_polymer *Allo_Fragment(int npol, int nback,int nangmax, FILE *flog)
{
	
	int ipol;
	struct s_polymer *f;
	//fprintf(stderr,"Allocando fragment lungo %d, nangmax %d\n",nback,nangmax);	

	f=(struct s_polymer *)malloc(npol*sizeof(struct s_polymer));
	if(!f)      Error("Cannot Allocate fragment");
	for(ipol=0;ipol<npol;ipol++)
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

void InvertTriang( double **Inv, double **Mat, int n )
{
 int i,k,a;

 for(i=0; i<n; i++)
   for(k=0; k<n; k++)
     Inv[i][k]=0.0;

 for(i=0; i<n; i++)
   Inv[i][i] = 1./Mat[i][i];

 for(i=0; i<n; i++)
  for(a=1; a<=i; a++)
  {
     for( k= i-a; k<i; k++)
       Inv[i][i-a] += (Mat[i][k] * Inv[k][i-a]);
       Inv[i][i-a] /= -Mat[i][i];
  }
return;
}




double Squared_n_Norma( double *vect, int dim )
{
 int i;
 double norm=0.0;

 for(i=0; i<dim; i++)
    norm += vect[i]*vect[i];

 return norm;
}

void TransposedMatOnVect( double **Mat, double *vect, double *risu, int n)
{
       int i,j;

      for(i=0; i<n; i++)
      {
            risu[i] = 0.0;
            for(j=0; j<n; j++)
                  risu[i] += Mat[j][i] * vect[j];
      }

      return;
}

int Gaussian_Angles(double *angles,int n)
{
        int i;
        double r1,r2;

      for(i=0;i<n;i++)
        {
                r1=(double)rand()/(double)RAND_MAX;
                r2=(double)rand()/(double)RAND_MAX;
                angles[i]=sqrt(-log(r1))*cos(2*PI*r2);
      }

      return 0;
}

void MatA( double **A, double **G, int dim)
{
      double a = PARAM_A;
      double b = PARAM_B;
      int i, j;

      for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
            {
                   
                   if(i==j)
                      A[i][i]= a*(1+ b*G[i][j])/2;
                  else
                        A[i][j]= a*b*G[i][j]/2;
            
      //fprintf(stderr,"A %d%d = %lf",i,j,A[i][j]);
              }
                    return;
}
      //
void Cholesky_2( double **L, double **A, int dim)
{
      int i,j,k;
      double sum;

      for( i=0; i<dim; i++)
            for( j=i; j<dim; j++)
            {
                  sum=A[i][j];
                  for(k=0;k<i;k++)
                        sum += -L[i][k]*L[j][k];
                  if(i == j)
                  {
                        if(sum<=0.0)
                              Error("Cholesky decomposition failed");
                        L[i][i]= sqrt(sum);
                  }
                  else
                        L[j][i]= sum/L[i][i];
            }

      return;
}
      

int CopyFragment(struct s_polymer *p,struct s_polymer *f,int iw,int nmul,int natom_fragment,int ip)
{
        int ok=1;
        int i;

      for(i=0;i<natom_fragment;i++)
        {
        CopyResidueCoordinates( (((p+ip)->back)+(iw-nmul+i)), (((f+ip)->back)+i));
        }

        return ok;
}


void CopyResidueCoordinates(struct s_back *from, struct s_back *to)
{
      (to)->nside=(from)->nside;
      (to->pos).x = (from->pos).x;
      (to->pos).y = (from->pos).y;
      (to->pos).z = (from->pos).z;
      (to)->move=(from)->move;
      strncpy(to->type,from->type,5);
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




      
int Compute_G(double **g,struct s_polymer *fragment,struct s_polymer *p,int iw,int nmul,int nang,struct s_mc_parms *parms)
{
        int ok=1;
	int natom_fragment=nmul+3;
        int i,j,l,i_ang;
	int ip=0;
      int k1;
      int k2;
      int k3;
      double deriv1[nang][3],deriv2[nang][3],deriv3[nang][3];
      double x1,y1,z1,x2,y2,z2,x3,y3,z3;
      double fact=(1./(LM_DELTAA));

      CopyFragment((p+ip),(fragment+ip),iw,nmul,natom_fragment,ip);

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
      for(i=0;i<natom_fragment;i++)
      {
            if( (((fragment+ip)->back)+i+1)->move==1)
            {
//		fprintf(stderr,"faccio il pivot %d \n",i);
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

                  CopyFragment((p+ip),(fragment+ip),iw,nmul,natom_fragment,ip);
                  i_ang++;
                  if(i_ang==(nang)) break;
            }
      }

      for(i=0;i<nang;i++)
            for(j=0;j<nang;j++)
                  for(l=0;l<3;l++)
      			{                  
			g[i][j]=(deriv1[i][l]*deriv1[j][l])+(deriv2[i][l]*deriv2[j][l])+(deriv3[i][l]*deriv3[j][l]);

			}

      return ok;
}



int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
	int ok=1,nang;
	int iw,ip;
	int i,m;
	double WN,WD;
	double deltaE;
		
	ip=0;
	iw=(nmul-1)+irand( (p+ip)->nback-nmul+1   ) ;
	nang=0;
	for(i=0;i<nmul;i++)
	{
		nang+=((((p+ip)->back)+iw-nmul+i)+1)->move;
	}

	deltaE = -GetEnergyMonomerRange(p,iw-nmul,iw,ip);

	if(iw==((p+ip)->nback-1)) //pivot forward OK
	{


		for(i=0;i<nang;i++)
			(fragment+ip)->d_ang[i]=parms->dw_mpivot*(0.5-frand());
		if (!parms->nodihpot)
                	for(i=0;i<nmul;i++)
				deltaE-=(((p+ip)->back)+iw-nmul+i)->e_dih;
		m=0;
		for(i=0;i<nmul;i++)
		{
			if( (((p+ip)->back)+iw+i-nmul+1)->move==1)
                	{
                     		ok*=PivotForward(p,iw-nmul+i,(fragment+ip)->d_ang[m],nmul-i,parms);
                     		m++;
                     		if(m==nang) break;
                	}

		}
		ok=AddSidechain(p,iw-nmul,iw,ip);
		deltaE += EnergyMonomerRange(p,pot,iw-nmul,iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
		if (!parms->nodihpot)
                	for(i=0;i<nmul;i++)
                        	deltaE+=EnergyDihedrals(p,pot,iw-nmul+i,ip,1);
		ok=Metropolis(deltaE,t,p->tables);

		if(ok==0)
                	UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
            	else
		{
                	UpdateMonomerRange(p,oldp,0,iw,ip,parms->shell);
                	p->etot+=deltaE;
                	parms->acc++;
		}


	} //end of forward
	
	
	else if(iw==nmul-1) //pivot backward OK
	{


		for(i=0;i<nang;i++)
                        (fragment+ip)->d_ang[i]=parms->dw_mpivot*(0.5-frand());
                if (!parms->nodihpot)
                        for(i=0;i<nmul;i++)
                                deltaE-=(((p+ip)->back)+i)->e_dih;
                m=0;
                for(i=0;i<nmul;i++)
                {
                        if( (((p+ip)->back)+iw+i-nmul+1)->move==1)
                        {
                                ok*=PivotBackward(p,iw-i,(fragment+ip)->d_ang[m],nmul-i-1,parms);
                                m++;
                                if(m==nang) break;
                        }

                }
                ok=AddSidechain(p,iw-nmul,iw,ip);
                deltaE += EnergyMonomerRange(p,pot,iw-nmul,iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
                if (!parms->nodihpot)
                        for(i=0;i<nmul;i++)
                                deltaE+=EnergyDihedrals(p,pot,iw-nmul+i,ip,1);
                ok=Metropolis(deltaE,t,p->tables);

                if(ok==0)
                        UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
                else
                {
                        UpdateMonomerRange(p,oldp,0,iw,ip,parms->shell);
                        p->etot+=deltaE;
                        parms->acc++;
                }
	
		
	} //end of backward


	else //pivot local
	{
		Compute_G((fragment+ip)->G,fragment,p,iw,nmul,nang,parms);	
		
		MatA((fragment+ip)->A,(fragment+ip)->G,nang);

		Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);

		Gaussian_Angles((fragment+ip)->g_ang,nang);
		InvertTriang((fragment+ip)->Y,(fragment+ip)->L,nang);
		TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang,nang);
		
		if(!parms->nodihpot)
			for(i=0;i<nmul;i++)
				deltaE -= (((p+ip)->back)+iw-nmul+i)->e_dih;

		if(!parms->noangpot)
		{
                	deltaE-=(((p+ip)->back)+iw-1)->e_ang;
                	deltaE-=(((p+ip)->back)+iw)->e_ang;
		}

		m=0;

		for(i=0;i<nmul;i++)
		{
			if( (((p+ip)->back)+iw+i-nmul+1)->move==1)
			{
				ok*=PivotForward(p,iw-nmul+i,(fragment+ip)->d_ang[m],nmul-i-2,parms);
				m++;
				if(m==nang)
				break;
                  	}
	
		}

		int out;
		double rc2=parms->r_cloose*parms->r_cloose;

	
	
		if( DAbs (Dist2( (((p+ip)->back)+iw-1)->pos,(((p+ip)->back)+iw)->pos ) -(((p+ip)->back)+iw-1)->d2_next < rc2 )  &&   DAbs( Angle( (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (((p+ip)->back)+iw+2)->pos, (p+ip)->tables, &out)  - (((p+ip)->back)+iw+1)->a_next) < 1.0)        
		{
			





			ok=AddSidechain(p,iw-nmul,iw,ip);
			deltaE += EnergyMonomerRange(p,pot,iw-nmul,iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
			if (!parms->nodihpot)
        	        	for(i=0;i<nmul;i++)
                              		deltaE+=EnergyDihedrals(p,pot,iw-nmul+i,ip,1);
                  
                  	if(!parms->noangpot)
                  	{
                        	deltaE+=EnergyAngles(p,pot,iw-1,ip,1);
                        	deltaE+=EnergyAngles(p,pot,iw,ip,1);
                  	}

                  	WD=1.0;
                  	for(i=0;i<nang;i++)
                        	WD=WD*(fragment+ip)->L[i][i];
                  	WD=WD*exp(-Squared_n_Norma((fragment+ip)->g_ang,nang));
			

                  	Compute_G((fragment+ip)->G,fragment,p,iw,nmul,nang,parms);
			       
 	          	MatA((fragment+ip)->A,(fragment+ip)->G,nang);
 		
			Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
                  	TransposedMatOnVect((fragment+ip)->L,(fragment+ip)->d_ang,(fragment+ip)->g_ang,nang);

		  	WN=1.0;
                  	for(i=0;i<nang;i++)
                       		WN=WN*(fragment+ip)->L[i][i];
                  	WN=WN*exp(-Squared_n_Norma((fragment+ip)->g_ang,nang));

                  	ok=B_Metropolis(deltaE,t,WN,WD,p->tables);

			
			if(ok==0) //move rejected
                  	{
                        	UpdateMonomerRange(oldp,p,iw-nmul,iw,ip,parms->shell);
                   	}
                  	else	//move accepted
                  	{
                        UpdateMonomerRange(p,oldp,iw-nmul,iw,ip,parms->shell);
                        p->etot+=deltaE;
                        parms->acc++;
                  	}

			
			
			
		}//end of loose condition

		else //if not loose pivot
        	{


 	        	UpdateMonomerRange(oldp,p,iw-nmul,iw,ip,parms->shell);
            	}



	}//end of local pivot


	parms->mov++;

	return ok;

}	

