#include "montegrappa.h"


struct s_polymer *AlloFragment(int npol, int nmul_local, FILE *flog)
{
      int ipol,natom_fragment;
      struct s_polymer *f;

      natom_fragment=(nmul_local+1)*3; 
      f=(struct s_polymer *)malloc(npol*sizeof(struct s_polymer));
      if(!f)      Error("Cannot Allocate fragment");
      for(ipol=0;ipol<npol;ipol++)
      {
            (f+ipol)->nback=natom_fragment;
            (f+ipol)->back=(struct s_back *)malloc(natom_fragment*sizeof(struct s_back));
            if(!(f+ipol)->back)     Error("\tAlloFragment(): Cannot Allocate backbone");     
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


int print_square_matrix(double **m,int n)
{
      int i,j;
      for(i=0;i<n;i++)
      {
            //fprintf(stderr,"\n");
            for(j=0;j<n;j++)
              //    fprintf(stderr,"(%d,%d)=%lf\t",i,j,m[i][j]);
                    fprintf(stderr,"%lf\t",m[i][j]);
      }
      
      return 0;
}//end of print

int print_vector(double *v,int n)
{
      int i;
      
      for(i=0;i<n;i++)
      {
           fprintf(stderr,"\n%lf\t\n",v[i]);
      }
      return 0;
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
            }
      
      return;
}

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



//copio un frammento lungo n a partire da k della catena ip
int CopyFragment(struct s_polymer *p,struct s_polymer *f,int offset,int n,int ip)
{
	int ok=1;
	int i;

      for(i=0;i<n;i++)
	{
      	CopyResidueCoordinates( (((p+ip)->back)+(offset+i)-(NMUL*3)), (((f+ip)->back)+i));
	}

	return ok;
}



int ComputeG(double **g,struct s_polymer *fragment,struct s_polymer *p,int k,int n,int np,struct s_mc_parms *parms)
{
	int ok=1;
      int natom_back=3*NMUL;
      int natom_fragment=(NMUL+NCHECK)*3;
	int i,j,l,i_ang;
	int ip=0;
      int k1;
      int k2;
      int k3;
      double deriv1[NMUL*2][3],deriv2[NMUL*2][3],deriv3[NMUL*2][3];
      double x1,y1,z1,x2,y2,z2,x3,y3,z3;
      double fact=(1./(LM_DELTAA));

      CopyFragment((p+ip),(fragment+ip),3*k,natom_fragment,ip);

      k1=natom_back-2;
      k2=natom_back-1;
      k3=natom_back;


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

                  deriv1[i_ang][0]=+(deriv1[i_ang][0]*fact);
                  deriv1[i_ang][1]=+(deriv1[i_ang][1]*fact);
                  deriv1[i_ang][2]=+(deriv1[i_ang][2]*fact);
                  deriv2[i_ang][0]=+(deriv2[i_ang][0]*fact);
                  deriv2[i_ang][1]=+(deriv2[i_ang][1]*fact);
                  deriv2[i_ang][2]=+(deriv2[i_ang][2]*fact);
                  deriv3[i_ang][0]=+(deriv3[i_ang][0]*fact);
                  deriv3[i_ang][1]=+(deriv3[i_ang][1]*fact);
                  deriv3[i_ang][2]=+(deriv3[i_ang][2]*fact);

                  CopyFragment(p,fragment,3*k,natom_fragment,ip);
                  i_ang++;
                  if(i_ang==(NMUL*2)) break;
            }   
      }

      for(i=0;i<i_ang;i++)
            for(j=0;j<i_ang;j++)
                  for(l=0;l<3;l++)
                        g[i][j]=(deriv1[i][l]*deriv1[j][l])+(deriv2[i][l]*deriv2[j][l])+(deriv3[i][l]*deriv3[j][l]);      

      return ok;
}



int MoveBiasedGaussian(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
	int ok=1;
	int iw,ip;
    	int i;
      double WN,WD;
	double deltaE;

	
      ip=0;
      iw=(NMUL)+irand( (((p+ip)->nback))/3-(NMUL) +1) ;
      deltaE = -GetEnergyMonomerRange(p,(iw*3)-(NMUL*3),(3*iw),ip); 

      //pivot forward
      if(iw==(((p+ip)->nback/3)))
      {
            int pippa;
            for(pippa=0;pippa<NMUL*2;pippa++){
                  (fragment+ip)->d_ang[pippa]=parms->dw_mpivot * (0.5 - frand());
            }
            if (!parms->nodihpot)
                  for (i=0;i<NMUL*3;i++) deltaE -= (((p+ip)->back)+(3*iw)-(NMUL*3)+i)->e_dih; 

            if(!parms->noangpot)
            {
                  deltaE-=(((p+ip)->back)+(3*iw)-(NMUL*3)-1)->e_ang;
                  deltaE-=(((p+ip)->back)+(3*iw)-(NMUL*3))->e_ang;
            }

       
            int m=0;
            for(i=0;i<(NMUL)*3;i++)
            {
                if( (((p+ip)->back)+(3*iw)+i-(NMUL*3)+1)->move==1)
                {
                     ok*=PivotForward(p,(3*iw)-(NMUL*3)+i,(fragment+ip)->d_ang[m],(NMUL*3)-i-2,parms);
                     m++;
                     if(m==NMUL*2) break;
                }

            }

            
            
            ok=AddSidechain(p,(3*iw)-(NMUL*3),(3*iw),ip);
            deltaE += EnergyMonomerRange(p,pot,(3*iw)-(NMUL*3),3*iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);     
            if (!parms->nodihpot)
                  for(i=0;i<NMUL*3;i++)
                        deltaE+=EnergyDihedrals(p,pot,(3*iw)-(NMUL*3)+i,ip,1);

            if(!parms->noangpot)
            {
                  deltaE+=EnergyAngles(p,pot,(3*iw)-(NMUL*3)-1,ip,1);
                  deltaE+=EnergyAngles(p,pot,(3*iw)-(NMUL*3),ip,1);
            }

            ok=Metropolis(deltaE,t,p->tables);

            if(ok==0)
                 UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
            else
            {
                 UpdateMonomerRange(p,oldp,0,(3*iw),ip,parms->shell);
                 p->etot+=deltaE;
                 parms->acc++;
            }

  
      }//end pivot forward
      

 
      //pivot backward     
      else if(iw==NMUL)
      {

            int pippa;
            for(pippa=0;pippa<NMUL*2;pippa++)
                  (fragment+ip)->d_ang[pippa]=parms->dw_mpivot * (0.5 - frand());
            
            if (!parms->nodihpot)
                  for (i=0;i<NMUL*3;i++)
                        deltaE -= (((p+ip)->back)+(3*iw)-(NMUL*3)+i)->e_dih;
            if(!parms->noangpot)
            {
                  deltaE-=(((p+ip)->back)+0)->e_ang;
                  deltaE-=(((p+ip)->back)+1)->e_ang;
            }


            int m=0;
            for(i=0;i<(NMUL)*3;i++)
            {
                  if( (((p+ip)->back)+(3*iw)+i-(NMUL*3)+1)->move==1)
                  {
                        ok*=PivotBackward(p,(3*iw)-i-1,(fragment+ip)->d_ang[m],(NMUL*3)-i-2,parms);
                        m++;
                        if(m==NMUL*2) break;    
                  }

            }

            ok=AddSidechain(p,(3*iw)-(NMUL*3),(3*iw),ip);

            deltaE += EnergyMonomerRange(p,pot,(3*iw)-(NMUL*3),3*iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);  
            if (!parms->nodihpot)
                  for(i=0;i<NMUL*3;i++)
                        deltaE+=EnergyDihedrals(p,pot,(3*iw)-(NMUL*3)+i,ip,1);
            if(!parms->noangpot)
            {
                  deltaE+=EnergyAngles(p,pot,0,ip,1);
                  deltaE+=EnergyAngles(p,pot,1,ip,1);
            }

            ok=Metropolis(deltaE,t,p->tables);

            if(ok==0)
                  UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
            else
            {
                 UpdateMonomerRange(p,oldp,0,(3*iw),ip,parms->shell);
                 p->etot+=deltaE;
                 parms->acc++;
            }

      }//end pivot backward



      //else, pivot locally
      else
      {
            ComputeG((fragment+ip)->G,fragment,p,iw,NMUL,1,parms);
            MatA((fragment+ip)->A,(fragment+ip)->G,NMUL*2);
            Cholesky_2((fragment+ip)->L,(fragment+ip)->A,NMUL*2);
            Gaussian_Angles((fragment+ip)->g_ang,NMUL*2);
            int m=0;
            InvertTriang((fragment+ip)->Y,(fragment+ip)->L,NMUL*2);
            TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang, NMUL*2);

            if (!parms->nodihpot)
                  for(i=0;i<NMUL*3;i++)
                        deltaE -= (((p+ip)->back)+(3*iw)-(NMUL*3)+i)->e_dih; 
            if(!parms->noangpot)
            {
                  deltaE-=(((p+ip)->back)+(3*iw)-(NMUL*3)-1)->e_ang;
                  deltaE-=(((p+ip)->back)+(3*iw)-(NMUL*3))->e_ang;
            }

            
            for(i=0;i<(NMUL)*3;i++)
            {
                  if( (((p+ip)->back)+(3*iw)+i-(NMUL*3)+1)->move==1)
                  {
                        ok*=PivotForward(p,(3*iw)-(NMUL*3)+i,(fragment+ip)->d_ang[m],(NMUL*3)-i-2,parms);
                        m++;
                        if(m==NMUL*2)
                              break;      
                  }

            }
            int out;
            double rc2=parms->r_cloose*parms->r_cloose;
      
            //loose pivot condition
            if( (Abs( Dist2( (((p+ip)->back)+(iw*3)-1)->pos,(((p+ip)->back)+((iw*3)))->pos ) - (((p+ip)->back)+(iw*3)-1)->d2_next) < rc2) && ( DAbs( Angle( (((p+ip)->back)+(iw*3)-2)->pos, (((p+ip)->back)+(iw*3)-1)->pos, (((p+ip)->back)+(iw*3))->pos, (p+ip)->tables, &out)  - (((p+ip)->back)+(iw*3)-1)->a_next) < 1.0))            
            {
                  ok=AddSidechain(p,(3*iw)-(NMUL*3),(3*iw),ip);
                  deltaE += EnergyMonomerRange(p,pot,(3*iw)-(NMUL*3),3*iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);        

                  if (!parms->nodihpot)
                  {                                                                               // new dihedral energy
                        for(i=0;i<NMUL*3;i++)
                              deltaE+=EnergyDihedrals(p,pot,(3*iw)-(NMUL*3)+i,ip,1);
                  }
                  if(!parms->noangpot)
                  {
                        deltaE+=EnergyAngles(p,pot,(3*iw)-(NMUL*3)-1,ip,1);
                        deltaE+=EnergyAngles(p,pot,(3*iw)-(NMUL*3),ip,1);
                  }

                  WD=1.0;
                  for(i=0;i<NMUL*2;i++)
                        WD=WD*(fragment+ip)->L[i][i];
                  WD=WD*exp(-Squared_n_Norma((fragment+ip)->g_ang,NMUL*2));

                  ComputeG((fragment+ip)->G,fragment,p,iw,NMUL,1,parms);
                  MatA((fragment+ip)->A,(fragment+ip)->G,NMUL*2);
                  Cholesky_2((fragment+ip)->L,(fragment+ip)->A,NMUL*2);
                  TransposedMatOnVect((fragment+ip)->L,(fragment+ip)->d_ang,(fragment+ip)->g_ang,NMUL*2);
            
                  WN=1.0;
                  for(i=0;i<NMUL*2;i++)
                        WN=WN*(fragment+ip)->L[i][i];
                  WN=WN*exp(-Squared_n_Norma((fragment+ip)->g_ang,NMUL*2));
                 
	          ok=B_Metropolis(deltaE,t,WN,WD,p->tables);
		  
                  if(ok==0)
		  {

                        UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
                   }  
		  else
                  {
                        UpdateMonomerRange(p,oldp,0,(3*iw),ip,parms->shell);
                        p->etot+=deltaE;
                        parms->acc++;
                  }

            }//loose pivot condition

            else //if not loose pivot
            {
                  UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,parms->shell);
            }
 
      }//end of local pivot

      parms->mov++;

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
