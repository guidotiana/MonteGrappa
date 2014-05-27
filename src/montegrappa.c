/*
 * montegrappa.c
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include <unistd.h>

int main(int argc, char *argv[])
{
	struct s_polymer *polymer, *startp, *oldp,*fragment,*replica;
	struct s_mc_parms *parms;
	struct s_potential *u;
	char polname[50],potname[50],fname[50];
	int i,j,nat,ntypes,irun,ci,out,nrestr;

	FILE *ftrj=NULL,*fe=NULL,*fproc=NULL;
	time_t rawtime;
	double **ematold;

	if (argc<4)
	{
		fprintf(stderr,"\nusage: MonteGrappa polymer.pol potential.pot parameters.par\n");
		exit(1);
	}


	int nprocs;

	int my_rank=0;
	#ifdef ACTIVE_MPI

	MPI_Init(&argc,&argv);
	MPI_Status astatus;
	struct s_mpi_parms *mpiparms=malloc(sizeof(struct s_mpi_parms));

	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);	

	Create_parms_datatype(&(mpiparms->Parms_mpi));
        Create_back_datatype(&(mpiparms->Back_mpi));
        Create_side_datatype(&(mpiparms->Side_mpi));
        Create_rot_datatype(&(mpiparms->Rot_mpi));
        Create_pot_datatype(&(mpiparms->Pot_mpi));
		
	mpiparms->my_rank=my_rank;









	mpiparms->nprocs=nprocs;
	#endif
	
	if(my_rank==0)
	{
		Welcome(stderr);
		strcpy(polname,argv[1]);
		strcpy(potname,argv[2]);
		parms = ReadMcParms(argv[3]);

		// if kill -SIGUSR1 given, soft exit
		signal(SIGUSR1,SoftExit);

		#ifndef DEBUG
     			if (parms->debug>0) Error("To debug, the program must be compiled with #define DEBUG");
    		#else
			fprintf(stderr,"WARNING: executable compiled with DEBUG option\n");
		#endif

		#ifdef ACTIVE_MPI
            	
		fprintf(stderr,"~MPI: %d Procs\n",nprocs);
            	if(nprocs!= parms->ntemp)
            	{
                	Error("WARNING! ntemp must be equal to nprocs");
                  	return -1;
		}
            	#endif

	}
	#ifdef ACTIVE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	#ifdef ACTIVE_MPI
	parms=send_parms(mpiparms->my_rank,nprocs,mpiparms->Parms_mpi,parms,astatus);
	if (strcmp(parms->fnproc,""))
      	{
        	sprintf(fname,"%s_%d.log",parms->fnproc,my_rank);
        	fproc = fopen(fname,"w");
            	if (!fproc)
			Error("Cannot open log-process file for writing");
	}

      	if(parms->debug>1) 
		fprintf(fproc, "iproc=%d\tnrun=%d\tnstep=%d\tnpol=%d\ttemp=%f\tnstep_exchange=%d\n", my_rank, parms->nrun, parms->nstep, parms->npol, parms->T[mpiparms->my_rank],parms->nstep_exchange);
      	#else
      	fproc = stderr;
      	#endif



	parms->seed = Randomize(parms->seed);

	startp = AlloPolymer(parms->npol,NRESMAX,NSIDEMAX,NROTMAX,NATOMMAX,parms->shell,parms->nosidechains,parms->flog);
	polymer = AlloPolymer(parms->npol,NRESMAX,NSIDEMAX,NROTMAX,NATOMMAX,parms->shell,parms->nosidechains,parms->flog);
	oldp = AlloPolymer(parms->npol,NRESMAX,NSIDEMAX,NROTMAX,NATOMMAX,parms->shell,parms->nosidechains,parms->flog);
	replica = AlloPolymer(parms->npol,NRESMAX,NSIDEMAX,NROTMAX,NATOMMAX,parms->shell,parms->nosidechains,fproc);

      	if( NMUL !=0)
     	{
        	fprintf(stderr,"Fragment length: %d(+1)\n",parms->nmul_local);
		int nmul_local=NMUL;      
	        fragment=AlloFragment(parms->npol,nmul_local,parms->flog);
        	for(i=0;i<parms->npol;i++)
           	{
                	(fragment+i)->A=AlloDoubleMatrix(nmul_local*2,nmul_local*2);
                	(fragment+i)->G=AlloDoubleMatrix(nmul_local*2,nmul_local*2);
                	(fragment+i)->L=AlloDoubleMatrix(nmul_local*2,nmul_local*2);
                	(fragment+i)->Y=AlloDoubleMatrix(nmul_local*2,nmul_local*2);
                	(fragment+i)->g_ang=AlloDouble(nmul_local*2);      
                	(fragment+i)->d_ang=AlloDouble(nmul_local*2);
            	}       
	}
	

	startp->tables = InitTables(fproc);                         // keep tables only in (polymer+0)
	for (i=0;i<parms->npol;i++) (polymer+i)->tables = startp->tables; // tables of all polymers point to same address
	for (i=0;i<parms->npol;i++) (replica+i)->tables = startp->tables; // tables of all polymers point to same address
	for (i=1;i<parms->npol;i++) (startp+i)->tables = startp->tables;  // tables of all polymers point to same address

	// Read polymer stuff from file
	if(my_rank==0)
	{
	
	ReadPolymer(polname,startp,parms->flog,parms->npol, parms->debug, &nat, &ntypes);
	
	}

	#ifdef ACTIVE_MPI
	
        for (i=0;i<parms->npol;i++)
	{
        	send_struct(&((startp+i)->nback),mpiparms-> my_rank, nprocs, &nat, &ntypes, astatus);
        	startp = send_pol(mpiparms->my_rank, nprocs, (startp+i)->nback, mpiparms->Back_mpi, mpiparms->Side_mpi, mpiparms->Rot_mpi, startp, astatus, i, parms->nshell, parms->nosidechains);
      	}
	if(parms->debug>1) fprintf(fproc,"iproc=%d\tnat=%d\tntypes=%d\tnback=%d\n",mpiparms->my_rank, nat, ntypes, startp->nback);

      	#endif

	for (i=0;i<parms->npol;i++) if (!AddSidechain(startp,0,startp->nback,i)) Error("Cannot add sidechain in starting chain");

	// Read potential (needs the number of atoms from ReadPolymer)
	u = AlloPotential(nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb);
	if(my_rank==0)
	{
	
		ReadPotential(potname,u,parms,nat,ntypes);
	
	}
	#ifdef ACTIVE_MPI
		
	u = send_pot(nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb, my_rank, nprocs,mpiparms-> Pot_mpi, u, astatus);
//	 fprintf(fproc,"iproc=%d\trhc=%f\tkr2splice=%f\tkesplice=%f\n",my_rank, u->g_r0hard, u->kr2_splice, u->ke_splice);
//	fprintf(fproc,"iproc=%d\trhc=%f\tkr2splice=%f\tkesplice=%f\tdih_ram=%d\n",my_rank, u->g_r0hard, u->kr2_splice, u->ke_splice, u->dih_ram);
		
	//POTENZIALE SUI DIEDRI 
	if(u->dih_ram)
      	{
        	for(i=0;i<2;i++) for(j=0;j<2;j++) if(u->sigma[i][j]<=0) if(my_rank==0) Error("Sigma <=0 !!");
        	if(u->dih0[0][0]!= -57) if(my_rank==0) fprintf(stderr,"WARNING! phi0_a is usually equal to -57!!!!!!!!\n");
        	if(u->dih0[1][0]!= -129) if(my_rank==0) fprintf(stderr,"WARNING! phi0_b is usually equal to -129!!!!!!!!\n");
       		if(u->dih0[0][1]!= -47) if(my_rank==0) fprintf(stderr,"WARNING! psi0_a is usually equal to -47!!!!!!!!\n");
        	if(u->dih0[1][1]!= 124) if(my_rank==0) fprintf(stderr,"WARNING! psi0_b is usually equal to 124!!!!!!!!\n");
     	}
      	
	#endif

	#ifdef OPTIMIZEPOT
	

	char fname_op[50];
	FILE *fp;
	#ifdef ACTIVE_MPI
	sprintf(fname_op,"op_log_proc%d.pdb",my_rank);
	#else
	sprintf(fname_op,"op_log");
	#endif
	fp = fopen(fname_op,"w");
	startp->op = InitializeOptimizePot(parms,ntypes,u,fproc,my_rank,fp,&nrestr);	
	ematold = AlloDoubleMatrix(ntypes,ntypes);

	#endif

	CheckOverlaps(startp,u,parms,parms->npol,parms->nosidechains,1,fproc);

	// calculate initial distances, angles and dihedrals to constrain loose moves
	if (parms->debug>1) fprintf(stderr,"Initial distances, angles and dihedrals: ");
	for (ci=0;ci<parms->npol;ci++)
		for (i=0;i<(startp+ci)->nback-1;i++)
		{
			(((startp+ci)->back)+i)->d2_next = Dist2( (((startp+ci)->back)+i)->pos, (((startp+ci)->back)+i+1)->pos );
			if (i>0) 
			{
				(((startp+ci)->back)+i)->a_next = Angle( (((startp+ci)->back)+i-1)->pos , (((startp+ci)->back)+i)->pos,
					(((startp+ci)->back)+i+1)->pos, startp->tables, &out);
//                        fprintf(stderr,"\n calcolo a_next = %lf",(((startp+ci)->back)+i)->a_next );
				if (!out) Error("Zero angle in initial structure");
			}
			if (i>1)
			{ 
				(((startp+ci)->back)+i)->d_next = Dihedral( (((startp+ci)->back)+i-2)->pos, (((startp+ci)->back)+i-1)->pos ,
					(((startp+ci)->back)+i)->pos, (((startp+ci)->back)+i+1)->pos, startp->tables, &out);
				if (!out) Error("Flat angle in initial structure");
			}
			if (parms->debug>1) fprintf(stderr,"d2(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->d2_next);
			if (parms->debug>1 && i>0) fprintf(stderr,"ang(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->a_next);
			if (parms->debug>1 && i>1) fprintf(stderr,"dih(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->d_next);
		}
	if (parms->debug>1) fprintf(stderr,"\n");

	fprintf(stderr,"Starting simulation...\n");

	// do Monte Carlo
	for (irun=0;irun<parms->nrun;irun++)
	{
		if (irun==0 || parms->always_restart)
		{
			CopyAllPolymers(startp,polymer,parms->npol,parms->nosidechains,parms->nosidechains);
			CopyAllPolymers(startp,replica,parms->npol,parms->nosidechains,parms->nosidechains);
		}
		polymer->etot = TotalEnergy(polymer,u,parms,parms->npol,1,parms->nosidechains,parms->debug,my_rank);
		replica->etot = TotalEnergy(replica,u,parms,parms->npol,1,parms->nosidechains,parms->debug,my_rank);
		if(my_rank==0)
		{
                	fprintf(stderr,"Energia=%lf\tstep=%d\n",polymer->etot,irun);
                	fprintf(stderr,"Energia_rep=%lf\tstep=%d\n",replica->etot,irun);
            	}

		if (parms->shell==1)	UpdateShell(polymer,parms);
		if (parms->nrun>1) fprintf(stderr,"NRUN = %d\n",irun);
		fprintf(stderr,"Proc %d: Initial energy = %lf\n",my_rank,polymer->etot);

		// open trajectory file
                  #ifdef ACTIVE_MPI
                  if (parms->nrun>1) sprintf(fname,"%s_run%d_proc%d.pdb",parms->fntrj,irun,my_rank);
                  else sprintf(fname,"%s_proc%d.pdb",parms->fntrj,my_rank);
                  #else
                  if (parms->nrun>1) sprintf(fname,"%s_run%d.pdb",parms->fntrj,irun);
                  else sprintf(fname,"%s.pdb",parms->fntrj);
                  #endif
                  ftrj = fopen(fname,"w");
                  if (!ftrj) Error("Cannot open trajectory file for writing");
                  sprintf(polymer->title,"initial\tE=%lf",polymer->etot);
                  PrintPDBStream(polymer,parms->npol,ftrj);

		

	// open energy file
		if (strcmp(parms->fne,""))
		{
			#ifdef ACTIVE_MPI

	        	if (parms->nrun>1)
				sprintf(fname,"%s_run%d_proc%d.ene",parms->fne,irun,my_rank);
                  	else
				sprintf(fname,"%s_proc%d.ene",parms->fne,my_rank);
                  	#else
                  	if (parms->nrun>1) 
				sprintf(fname,"%s_run%d.ene",parms->fne,irun);
                  	else
				sprintf(fname,"%s.ene",parms->fne);

                  	#endif

			fe = fopen(fname,"w");
			if (!fe) 
				Error("Cannot open energy file for writing");
			PrintEnergies(fe,parms->npol,0,polymer);
		}

		

		// MAKE MONTE CARLO
		#ifdef ACTIVE_MPI		
		Do_MC(polymer,fragment,replica,startp,u,parms,ftrj,fe,oldp,fproc,irun,mpiparms);
		#else
		Do_MC(polymer,fragment,replica,startp,u,parms,ftrj,fe,oldp,fproc,irun,0);
		#endif

		 #ifdef OPTIMIZEPOT
	         	if (strcmp(parms->op_minim,"none"))
            		{
            			#ifdef ACTIVE_MPI       
            			
				if(my_rank==0) 
					OP_SamplePotential(polymer,parms,ntypes,u->e,irun,ematold);
            			u = send_pot(nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb, my_rank, nprocs, mpiparms->Pot_mpi, u, astatus);
            			MPI_Barrier(MPI_COMM_WORLD);

            			#else
                  		
				OP_SamplePotential(polymer,parms,ntypes,u->e,irun,ematold);
            
				#endif
				if(my_rank==0)
               			{
                			sprintf(fname,"newpot_%d.pot",irun);
                        		PrintPotential(u,fname,nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb);
                		}
            		}
            	#endif



		// open last conformation file
        	#ifdef ACTIVE_MPI
             	
		if (parms->nrun>1) 
			sprintf(fname,"%s_run%d_proc%d.pol",parms->flastp,irun,my_rank);
             	else
			sprintf(fname,"%s_proc%d.pol",parms->flastp,my_rank);
            	
		#else
             
		if (parms->nrun>1)
			sprintf(fname,"%s_run%d.pol",parms->flastp,irun);
             	else 
			sprintf(fname,"%s.pol",parms->flastp);
            	
		#endif
   
	        CheckOverlaps(polymer,u,parms,parms->npol,parms->nosidechains,1,fproc);
            	PrintPolymer(fname,polymer,parms->npol);
           	fprintf(fproc,"final: Etot=%lf\tEtot(true)=%lf\n",polymer->etot,TotalEnergy(polymer,u,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank));

		// close files
		if (strcmp(parms->fntrj,"")) fclose(ftrj);
		if (strcmp(parms->fne,"")) fclose(fe);

		fprintf(stderr,"End %d iteration\n",irun);
		#ifdef ACTIVE_MPI
            	if (my_rank==0) fprintf(stderr,"End %d iteration\n",irun);
            	#endif

	}

	#ifdef OPTIMIZEPOT
	fclose(fp);
        #endif


	free(polymer);
	free(oldp);
	free(ematold);
	free(startp);
	free(replica);
      	
	if(parms->nmul_local !=0 && parms->movetype[7]!=-1)
      {

      for(i=0;i<parms->npol;i++)
           free( (fragment+i)->back );
            
            for(j=0;j<parms->nmul_local*2;j++)
            {
                  free( (fragment+i)->A[j] );
                  free( (fragment+i)->G[j] );
                  free( (fragment+i)->L[j] );
                  free( (fragment+i)->Y[j] );
           }

      free( (fragment+i)->A );
      free( (fragment+i)->G );
      free( (fragment+i)->L );
      free( (fragment+i)->Y );
      free( (fragment+i)->d_ang );
      free( (fragment+i)->g_ang );


      free(fragment);
	

      }

	free(parms);

	
      #ifdef ACTIVE_MPI
      time(&rawtime);
      fprintf(fproc,"Run finished at %s\n\n",asctime(localtime(&rawtime)));
      MPI_Type_free(&mpiparms->Parms_mpi);
      MPI_Type_free(&mpiparms->Pot_mpi);
      MPI_Type_free(&mpiparms->Rot_mpi);
      MPI_Type_free(&mpiparms->Side_mpi);
      MPI_Type_free(&mpiparms->Back_mpi);


      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      #endif


	return 0;
}

void Welcome(FILE *fp)
{
/*
      fprintf(fp,"\n\n***************************************\n");
      fprintf(fp,"*        MonteGrappa v. %d.%d           *\n",NVER,NSUBVER);
      fprintf(fp,"***************************************\n");
      fprintf(fp,"  G. Tiana, 2010\n");
      fprintf(fp,"pid = %d\n",getpid());
*/

/*
fprintf(fp,"\n\n");
fprintf(fp,"  __  __             _              \n");
fprintf(fp," |  \\/  | ___  _ __ | |_ ___        \n");
fprintf(fp," | |\\/| |/ _ \\| '_ \\| __/ _ \\       \n");
fprintf(fp," | |  | | (_) | | | | ||  __/       \n");
fprintf(fp," |_|  |_|\\___/|_| |_|\\__\\___|       \n");
fprintf(fp,"                                    \n");
fprintf(fp,"   ____                             \n");
fprintf(fp,"  / ___|_ __ __ _ _ __  _ __   __ _ \n");
fprintf(fp," | |  _| '__/ _` | ' _\\| '_ \\ / _` |\n");
fprintf(fp," | |_| | | | (_| | |_) | |_) | (_| |\n");
fprintf(fp,"  \\____|_|  \\__,_| .__/| .__/ \\__,_|\n");
fprintf(fp,"                 |_|   |_|          \n");
fprintf(fp,"\nG. Tiana, 2010\n");
fprintf(fp,"pid = %d\n",getpid());
fprintf(fp,"\n\n");
*/



fprintf(fp,"\n\n");
fprintf(fp,"             ,/k.\n");
fprintf(fp,"            /  ih,\t\t*MonteGrappa*\n");     
fprintf(fp,"       ,-' ,  `:7b \t\t\tv1.0\n");
fprintf(fp,"     _.-/   '  /b.`.4p,\n");
fprintf(fp,"  --   ,    ,-' ^6x, `.'^=._\n");
fprintf(fp,"\n");
fprintf(fp,"\nG. Tiana, 2010\n");
fprintf(fp,"pid = %d\n",getpid());
fprintf(fp,"\n\n");



}

