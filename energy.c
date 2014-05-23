/*
 * energy.c
 *
 *  Created on: Oct 11, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include "grappino.h"

/*****************************************************************************
 Go two-body potential (in this case ntypes=natoms)
 *****************************************************************************/
void Go_Pairs(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int nchains, int nat, double **cm)
{
	int i,j,k=0;
	double r;
	FILE *fout;

	if (parms->go_noise) 
	{
		fprintf(stderr,"Adding noise (kind %d, sigma=%lf) to the Go potential\n",parms->go_noise,parms->go_noise_sigma);
		srand(time(0));	
	}
	if (parms->debug>1)  fprintf(stderr,"\n\nGo contacts:\n");

	for (i=0;i<nat;i++)
		for (j=i;j<nat;j++)
		{
			if (cm[i][j]>0)									// if it is a native contact
			{
				e[i][j] = parms->goe;
				if (parms->go_noise) e[i][j] += gauss(0.,parms->go_noise_sigma);	// if you want to add noise to the go potential
				e[j][i] = e[i][j];							// make the matrix symmetric
				if (parms->use_nativedist) r = cm[i][j] * cm[i][j] * parms->k_native_r * parms->k_native_r;
				else r = parms->rnat * parms->rnat;

				if (r2[i][j]>r) r = r2[i][j];					// if already set, use the maximum between the two
			
				r2[i][j] = r;
				r2[j][i] = r2[i][j];

				if (parms->use_nativedist) r = cm[i][j] * cm[i][j] * parms->k_native_hc * parms->k_native_hc;
				else r = parms->rhard * parms->rhard;

				if (r02[i][j]>0 && r02[i][j]<r) r=r02[i][j];
				r02[i][j] = r;
				r02[j][i] = r02[i][j];
				if (parms->debug>1)  fprintf(stderr,"%d-%d\t",i,j);
				k++;
			}

			if (parms->go_noise==2)								// if you want to add non-go noise
			{
				e[i][j] = gauss(0.,parms->go_noise_sigma);
				e[j][i] = e[i][j];
				r2[i][j] =  parms->rnat * parms->rnat;
				r2[j][i] = r2[i][j];
				r02[i][j] = parms->rhard * parms->rhard;
				r02[j][i] = r02[i][j];
			}
		}
	if (parms->debug>1)  fprintf(stderr,"\n");


	// if required, print native contacts to contactfile
	if (strcmp(parms->cntfile,""))
	{
		fout = fopen(parms->cntfile,"w");
		if (!fout) Error("Cannot open file to write native contacts");
		for (i=0;i<nat;i++)
			for (j=0;j<nat;j++)
				if (cm[i][j]>0) fprintf(fout,"%4d %4d\n",i,j);


		fclose (fout);
	}
	if (parms->debug>0) fprintf(stderr,"Found %d Go contacts\n",k);

}

void Go_Dihedrals(struct s_parms *parms, struct s_polymer *p, int nc, double *dih01, double *dih03, struct s_potential *u)
{
	int i,j,ia,out;

	for (i=0;i<nc;i++)
		for (j=2;j<(p+i)->nback-1;j++)
		{
			ia =  (((p+i)->back)+j)->ia;
			u->dih01[ia] = Dihedral( (((p+i)->back)+j-2)->pos, (((p+i)->back)+j-1)->pos, 
						(((p+i)->back)+j)->pos, (((p+i)->back)+j+1)->pos, p->tables,&out);
			u->e_dih1[ia] = parms->e_dih1;
			u->dih03[ia] = u->dih01[ia];
			u->e_dih3[ia] = parms->e_dih3;
		}
}

void Go_Angles(struct s_parms *parms, struct s_polymer *p, int nc, double *ang, struct s_potential *u)
{
	int i,j,ia,out;

	for (i=0;i<nc;i++)
		for (j=1;j<(p+i)->nback-1;j++)
		{
			ia =  (((p+i)->back)+j)->ia;
			u->ang0[ia] = Angle( (((p+i)->back)+j-1)->pos, (((p+i)->back)+j)->pos, (((p+i)->back)+j+1)->pos, p->tables,&out );
			u->e_ang[ia] = parms->e_ang;
		}
}


/*****************************************************************************
 Contact map between atoms sticking out of different backbone atoms
 (nat x nat)
 *****************************************************************************/
double **ContactMap(struct s_parms *parms, struct s_polymer *p, int nchains, int nat,int debug)
{
	double **cm;
	int ib1,ib2,is1,is2,c1,c2;
	double d;

	cm = AlloDoubleMatrix(nat,nat);
	if (debug>2) fprintf(stderr,"Contact map:\n");

	// loop on backbone atoms
	for (c1=0;c1<nchains;c1++)
		for (c2=c1;c2<nchains;c2++)
			for (ib1=0;ib1<(p+c1)->nback;ib1++)
				for (ib2=0;ib2<(p+c2)->nback;ib2++)
					if (c1!=c2 || ib1<=ib2-parms->imin)
					{
						// backbone-backbone interactions
						d = Dist( (((p+c1)->back)+ib1)->pos, (((p+c2)->back)+ib2)->pos );
						if ( d < parms->rnat )
						{
							cm[(((p+c1)->back)+ib1)->itype][(((p+c2)->back)+ib2)->itype] = d;
							cm[(((p+c2)->back)+ib2)->itype][(((p+c1)->back)+ib1)->itype] = d;
							if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((p+c1)->back)+ib1)->ia,(((p+c1)->back)+ib1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((p+c2)->back)+ib2)->ia,(((p+c2)->back)+ib2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);
						}

						// sidechain-backbone interactions
						for (is1=0;is1<(((p+c1)->back)+ib1)->nside;is1++)
						{
							d = Dist( (((((p+c1)->back)+ib1)->side)+is1)->pos, (((p+c2)->back)+ib2)->pos );
							if ( d < parms->rnat )
							{
								cm[(((((p+c1)->back)+ib1)->side)+is1)->itype][(((p+c2)->back)+ib2)->itype] = d;
								cm[(((p+c2)->back)+ib2)->itype][(((((p+c1)->back)+ib1)->side)+is1)->itype] = d;
								if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((((p+c1)->back)+ib1)->side)+is1)->ia,(((((p+c1)->back)+ib1)->side)+is1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((p+c2)->back)+ib2)->ia,(((p+c2)->back)+ib2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

							}
						}
						for (is2=0;is2<(((p+c2)->back)+ib2)->nside;is2++)
						{
							d = Dist( (((((p+c2)->back)+ib2)->side)+is2)->pos, (((p+c1)->back)+ib1)->pos );
							if ( d < parms->rnat )
							{
								cm[(((((p+c2)->back)+ib2)->side)+is2)->itype][(((p+c1)->back)+ib1)->itype] = d;
								cm[(((p+c1)->back)+ib1)->itype][(((((p+c2)->back)+ib2)->side)+is2)->itype] = d;
								if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((p+c1)->back)+ib1)->ia,(((p+c1)->back)+ib1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((((p+c2)->back)+ib2)->side)+is2)->ia,(((((p+c2)->back)+ib2)->side)+is2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

							}
						}
						// sidechain-sidechain interactions
						for (is1=0;is1<(((p+c1)->back)+ib1)->nside;is1++)
							for (is2=0;is2<(((p+c2)->back)+ib2)->nside;is2++)
							{
								d = Dist( (((((p+c1)->back)+ib1)->side)+is1)->pos, (((((p+c2)->back)+ib2)->side)+is2)->pos );
								if ( d < parms->rnat )
								{
									cm[(((((p+c1)->back)+ib1)->side)+is1)->itype][(((((p+c2)->back)+ib2)->side)+is2)->itype] = d;
									cm[(((((p+c2)->back)+ib2)->side)+is2)->itype][(((((p+c1)->back)+ib1)->side)+is1)->itype] = d;
									if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((((p+c1)->back)+ib1)->side)+is1)->ia,(((((p+c1)->back)+ib1)->side)+is1)->type,
										(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((((p+c2)->back)+ib2)->side)+is2)->ia,(((((p+c2)->back)+ib2)->side)+is2)->type,
										(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

								}
							}
					}

			return cm;
}

/*****************************************************************************
 Assign as atom types the atom id (and ntypes=natoms), as it is in the Go model
 *****************************************************************************/
int SetGoTypes(struct s_polymer *p, int nchains, int nat)
{
	int ic,i,j;

	fprintf(stderr,"Set Go types\n");

	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			(((p+ic)->back)+i)->itype = (((p+ic)->back)+i)->ia;
			for (j=0;j<(((p+ic)->back)+i)->nside;j++) (((((p+ic)->back)+i)->side)+j)->itype = (((((p+ic)->back)+i)->side)+j)->ia;
		}

	return nat;
}

/*****************************************************************************
 Read atom types from a library file
 *****************************************************************************/
int ReadTypes(struct s_polymer *p, int nchains, int nat, char *nfile)
{
	int i,ic,j,k,it=0,type[NTYPEMAX],nt=0;
	char aux[500],atom[NTYPEMAX][5],aa[NTYPEMAX][5];
	FILE *fp;

	// read from file
	fprintf(stderr,"Reading atom types from file %s\n",nfile);

	fp = fopen(nfile,"r");
	if (!fp) Error("Cannot open file for reading atom types, check atomtypes directive");

	while(fgets(aux,500,fp)!=NULL)
	{
		if ( sscanf(aux,"%d %s %s",&(type[it]),atom[it],aa[it]) == 3) it++;
		if (it>=NTYPEMAX) Error("NTYPEMAX too small in ReadTypes");
		if (it>nt) nt=it;
	}	

	fprintf(stderr,"Read %d atomtypes\n",it);
	fclose(fp);

	// assign types
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			// backbone
			(((p+ic)->back)+i)->itype = -1;
			for (k=0;k<it;k++)
				if (  ( !strcmp((((p+ic)->back)+i)->type, atom[k]) || !strcmp( "*", atom[k]) )
					&&  ( !strcmp( (((p+ic)->back)+i)->aa, aa[k]) || !strcmp( "*", aa[k])) )
						(((p+ic)->back)+i)->itype = type[k];
			if ( (((p+ic)->back)+i)->itype == -1)
			{
				fprintf(stderr,"Cannot assign atomtype to chain %d, back=%d (aa=%s atom=%s)\n",ic,i,
					(((p+ic)->back)+i)->aa,(((p+ic)->back)+i)->type);
				Error("Cannot assign atom type");
			}

			// sidechain
			for (j=0;j<(((p+ic)->back)+i)->nside;j++)
			{
				(((((p+ic)->back)+i)->side)+j)->itype = -1;
				for (k=0;k<it;k++)
					if ( ( !strcmp((((((p+ic)->back)+i)->side)+j)->type, atom[k]) || !strcmp( "*", atom[k]) )
						&& ( !strcmp( (((p+ic)->back)+i)->aa, aa[k]) || !strcmp( "*", aa[k])) )
							(((((p+ic)->back)+i)->side)+j)->itype = type[k];


				if ( (((((p+ic)->back)+i)->side)+j)->itype == -1)
				{
					fprintf(stderr,"Cannot assign atomtype to chain %d, back=%d (aa=%s atom=%s)\n",ic,i,
						(((p+ic)->back)+i)->aa,(((((p+ic)->back)+i)->side)+j)->type);
					Error("Cannot assign atom type");
				}
			}
		}

	return nt;
}


/*****************************************************************************
 Add a cys-cys interaction
 *****************************************************************************/
void DisulphideBonds(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int nchains, int nat)
{
	int i,j,ic,jc,ki,kj,it1,it2,cnt=0;

	fprintf(stderr,"Adding disulphide bonds...\n");

	for (ic=0;ic<nchains;ic++)
		for (jc=0;jc<nchains;jc++)
			for (i=0;i<(p+ic)->nback;i++)
				for (j=0;j<(p+jc)->nback;j++)
					if ( !strcmp( (((p+ic)->back)+i)->aa , "CYS" ) )
						if ( !strcmp( (((p+jc)->back)+j)->aa , "CYS" ) )
							for (ki=0;ki<(((p+ic)->back)+i)->nside;ki++)
								for (kj=0;kj<(((p+jc)->back)+j)->nside;kj++)
									if ( !strncmp( (((((p+ic)->back)+i)->side)+ki)->type , "S" ,1) )
										if ( !strncmp( (((((p+jc)->back)+j)->side)+kj)->type , "S" ,1) )
											if (ic!=jc || i!=j)
												{
													it1 = (((((p+ic)->back)+i)->side)+ki)->itype;
													it2 = (((((p+jc)->back)+j)->side)+kj)->itype;
													e[it1][it2] = parms->cys;
													e[it2][it1] = e[it1][it2];
													r2[it1][it2] = parms->cysr * parms->cysr;
													r2[it2][it1] = r2[it1][it2];
													if (r02[it1][it2]<EPSILON)
													{
														r02[it1][it2] = parms->rhard;
														r02[it2][it1] = r02[it1][it2];
													}
													cnt++;
												}
	fprintf(stderr,"Added %d disulphide bond interactions.\n",cnt);
}

void HydrophobicInteraction(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int natypes)
{
	int i,j,cnt=0;

	fprintf(stderr,"Adding hydrophobic interaction...\n");

	for (i=0;i<natypes;i++)
		for (j=0;j<natypes;j++)
			if (parms->hydro_at[i]==1 && parms->hydro_at[j]==1)
			{
				e[i][j] = parms->hydro_e;
				e[j][i] = e[i][j];
				r2[i][j] = parms->hydro_r * parms->hydro_r;
				r2[j][i] = r2[i][j];
				if (r02[i][j]<EPSILON)
				{
					r02[i][j] = parms->rhard;
					r02[j][i] = r02[i][j];
				}
				cnt++;
			}
}

/**********************************************
 Create OP file
 **********************************************/
void PrintOpGoFile(char *nfile, double **cm, struct s_polymer *p, int nchains, char *kind)
{
	int ic,jc,i,j,k,r;
	FILE *fp;

	fprintf(stderr,"Write OP file %s for potential optimization (kind=%s)\n",nfile,kind);

	fp = fopen(nfile,"w");
	if (!fp) Error("Cannot open op file for writing");

	// distances between native CA (numbers are ia)
	if (!strcmp(kind,"GO_DIST_CA"))
	{
		for (ic=0;ic<nchains;ic++)
			for (i=0;i<(p+ic)->nback;i++)	
				for (jc=0;jc<nchains;jc++)
					for (j=0;j<(p+jc)->nback;j++)
						if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") )
							fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((p+jc)->back)+j)->ia,
								Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos ) );

	}

	// contacts between backbones (numbers are iback, only one chain)
	if (!strcmp(kind,"GO_DIST_CONT"))
	{
		for (i=0;i<(p+0)->nback;i++)	
			for (j=i+4;j<(p+0)->nback;j++)
				for (k=0;k<(((p+0)->back)+i)->ncontacts;k++)
					if ( *(((((p+0)->back)+i)->contacts)+k) == j )
						fprintf(fp,"%d\t%d\t\t0\t1.\t0.2\n",i,j);
	}

	
		// distances between some native CA (with different sigma)
	if (!strcmp(kind,"GO_DIST_CA_2"))
	{
		int nrest = 0;
		double **restrains;
		int maxrest = ((p->nback)*(p->nback)*10);
		restrains = AlloDoubleMatrix(4,maxrest);
		
		for (ic=0;ic<nchains;ic++)
			for (i=0;i<(p+ic)->nback;i++)	
				for (jc=0;jc<nchains;jc++)
					for (j=i;j<(p+jc)->nback;j++)
						if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") )
						{
							double d = Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos );
							//fprintf(stderr, "dist = %lf\n", d);
							if(Abs(i-j)>=20 && d>0)
							{
								restrains[0][nrest]= (double) (((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+ic)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 1; 
								nrest++;
							}
							/*if(Abs(i-j)>=20){

								restrains[0][nrest]= (double)(((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+ic)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 1.5; 
								nrest++;
							}
							if(d>10 && Abs(i-j)>4 && Abs(i-j)<20){
								restrains[0][nrest]= (double)(((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+ic)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 3; 
								nrest++;
							}*/
							
						if(nrest>=maxrest) fprintf(stderr, "number of restrains too large\n");
						}
		fprintf(fp,"ndata %d\n",nrest);
		for (r=0;r<nrest;r++)
			fprintf(fp,"%g\t%g\t\t1\t%lf\t%f\n", restrains[0][r], restrains[1][r], restrains[2][r], restrains[3][r] );
		free(restrains);
	}


	//distances between some native CA, with restraints chosen from vmd
	if (!strcmp(kind,"GO_DIST_CA_3"))
	{
		int nrest = 0;
		double **restrains;
		int maxrest = ((p->nback)*(p->nback)*10);
		restrains = AlloDoubleMatrix(4,maxrest);
	
	
		for (ic=0;ic<nchains;ic++)
			for (i=0;i<(p+ic)->nback;i++)	
				for (jc=0;jc<nchains;jc++)
					for (j=i;j<(p+jc)->nback;j++)
						if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") )
						{
							double d = Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos );
							//fprintf(stderr, "dist = %lf\n", d);
							if(Abs(i-j)>=4 && Abs(i-j)<=5 && d>0)
							{
								restrains[0][nrest]= (double) (((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+jc)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 1; 
								nrest++;
							}
							if(Abs(i-j)>5 && d<=10 && d>0)
							{
								restrains[0][nrest]= (double)(((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+jc)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 1; 
								nrest++;
							}	
							if((i==2 && j==14 && d>10)||(i==2 && j==35 && d>10)||(i==2 && j==29 && d>10)||(i==0 && j==25 && d>10)||(i==14 && j==25 && d>10)||(i==14 && j==35 && d>10))
							{	
								restrains[0][nrest]= (double)(((p+ic)->back)+i)->ia;
								restrains[1][nrest]= (double)(((p+jc)->back)+j)->ia;
								restrains[2][nrest]= d;
								restrains[3][nrest]= 1; 
								nrest++;
							}
						if(nrest>=maxrest) fprintf(stderr, "number of restrains too large\n");
						}

		fprintf(fp,"ndata %d\n",nrest);
		for (r=0;r<nrest;r++)
			fprintf(fp,"%g\t%g\t\t1\t%lf\t%f\n", restrains[0][r], restrains[1][r], restrains[2][r], restrains[3][r] );
		free(restrains);
	}

	if(!strcmp(kind,"GO_DIST_ALLATOM"))		//ATTT! non considero mai l'ossigeno legato a C
	{
		int nrest = 0;
		double **restrains;
		int maxrest = ((p->nback)*(p->nback)*100);
		restrains = AlloDoubleMatrix(4,maxrest);
		int naa= p->nback/3;
		fprintf(stderr, "naa=%d\n",naa);
		int a,b;

		//distances between some backbone atoms
		for (i=0;i<naa;i++)	
			for (j=i;j<naa;j++)
				{
				a=3*i+rand()%3;
				b=3*j+rand()%3;
				double d = Dist( (((p)->back)+a)->pos, (((p)->back)+b)->pos );
				if(Abs(a-b)>8)
					{
					restrains[0][nrest]= (double)(((p)->back)+a)->ia;
					restrains[1][nrest]= (double)(((p)->back)+b)->ia;
					restrains[2][nrest]= d;
					restrains[3][nrest]= 0.4; 
					nrest++;
					}
				if(nrest>=maxrest) fprintf(stderr, "number of restrains too large\n");	
				}

		//distances between sidechain atoms
		for (i=0;i<naa;i++)	
			for (j=i+1;j<naa;j++){
				if(strcmp(((((p)->back)+1+3*i)->aa),"GLY") && strcmp(((((p)->back)+1+3*j)->aa),"GLY"))
				{
					a = rand()%((((p)->back)+1+3*i)->nside);
					b = rand()%((((p)->back)+1+3*j)->nside);
					double d = Dist( ((((p->back)+1+3*i)->side)+a)->pos, ((((p->back)+1+3*j)->side)+b)->pos );
					//if(Abs(i-j)<4||Abs(i-j)>15) per O_DIST_ALLATOM isata di solito
					if(Abs(i-j)<6||Abs(i-j)>20)
					{
						restrains[0][nrest]= (double) ((((p->back)+1+3*i)->side)+a)->ia;
						restrains[1][nrest]= (double) ((((p->back)+1+3*j)->side)+b)->ia;
						restrains[2][nrest]= d;
						restrains[3][nrest]= 0.6; 
						nrest++;
					}
				
				}
				if(nrest>=maxrest) fprintf(stderr, "number of restrains too large\n");	
			}

		fprintf(fp,"ndata %d\n",nrest);		
		for (r=0;r<nrest;r++)
			fprintf(fp,"%g\t%g\t\t1\t%lf\t%f\n", restrains[0][r], restrains[1][r], restrains[2][r], restrains[3][r] );
		free(restrains);

	}	

	fclose(fp);
	

}
