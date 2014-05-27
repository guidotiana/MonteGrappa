/*
 * stempering.c
 *
 *  Created on: Nov 24, 2009
 *      Author: guido
 */

#include "stempering.h"

// Stempering is carried out if there is the file _STEMPERING in the current dir.
//
// The file contains:
// method [stempering | adaptive]		whether to carry out standard or adaptive stempering
// nstep <int>							tries a temperature change every nstep
// nprint <int>							prints output every nprint*nstep steps
// preamble <int>						do some constant-T steps before anything else
// ntemp <int>							for standard stempering, how many temperatures to use,
//										for adaptive, how many temperatures to start with
// debug <int>							0=silent, 1=say what's doing, 2=give the details,
//										3=be very verbose, 4=every cough
// temperatures							followed by list of temperatures and maybe weights
// <double> [<double>]
// <double> [<double>]
// ...
// tfile <filename>						output file for temperatures
// thefile <filename>					output file for thermodynamics
// dosfile <filename>					output file for density of states
// dumbfile <filename>					output file for dumb thermodynamics
// histofile <filename>					output file for histograms
// nprntt <int>							print thermodynamics every nprntt steps
// nprndumb <int>						print dumb averages every nprndumb steps
//
// only for adaptive:
// ntmax <int>							maximum number of temperatures allowed
// nsadj <int>							adapts temperatures/weights every nsadj
// emin <double>						minimum energy for the histogram
// emax <double>						maximum energy for the histogram
// ebin <double>						binning of the histogram
// anneal [sigma|prob|none]				how to set lower temperature (i.e. sigma = proportionally to
//										energy stdev; prob = at fixed exchange probability)
// k_new <double>						the lowest temperature is set such that its energy lies
//										k_new times the standard deviation from the next-to-lowest
// lp_new <double>						the lowest temperature is set such that the estimated exchange
//										log-probability is p_new
// lpthresh <double>					log of minimum probablity to remove a temperature
//										set 0 to disable, set 9 to calculate automatically
// hthresh <double>						threshold on overall probability to keep a histogram is htresh/ntemp (default 0.10)
// binthresh <double>					discard bin in histogram if counts < binthresh
// restart								restart from RESTART_ST file
// t_only_add <double>					if lowest T is below this, do not try to remove replicas
// keepall								uses all past histograms to calculate thermodynamics
// nkeepold <int>						if keepall active, keep only histograms of last %d run
// sumoldhisto							1=sum together histograms with same temperature (only if keepall),
//										2=delete older histograms with the same temperature
// removet								0=never remove old histogram, 99=remove any number, else=max # of temperatures to remove
// tstop								stop the simulation when this temperature is reached
// tnorm								when this temperature is reached start normal stempering
// deltat								in mhistogram discard histograms which are closer than deltat in temperature
// ttarget						when attempts to go below this, choose exactly this and then stop
// pdbnf						file to print pdb output
 
/*************************************************************
 * Initialize simulated tempering.
 * 		define a struct st_stuff *p in the main code
 * 		call InitStempering before the main loop
 *************************************************************/
struct st_stuff *InitSTempering(void)
{
    int i,r;
	double input_temp[NTEMPMAX],input_g[NTEMPMAX],dumb;
	struct st_stuff *p;
	char aux[500],aux_gm[500];
	FILE *fin;

	// do simulated tempering if it finds the params file

	fin = fopen("STEMPERING","r");
	if (!fin) return (struct st_stuff *)NULL;

    	// Alocate structure

	p = (struct st_stuff *) calloc(1,sizeof(struct st_stuff));

	// Read params file

	strcpy(p->nftout,"temperatures.dat");
	strcpy(p->nfdos,"dos.dat");
	strcpy(p->nfthe,"");
	strcpy(p->nfdumb,"");
	strcpy(p->nfdumb2,"");
	strcpy(p->nfhisto,"");
	strcpy(p->anneal,"");
	strcpy(p->nfthe,"thermdodyn.dat");
	strcpy(p->nftout,"temperatures.dat");
	strcpy(p->nfdumb,"dumb.dat");
	strcpy(p->pdbnf,"harvest.pdb");

	p->nprint = 1;
	p->k = 1.;
	p->p_new = 1.;		// >0 means not using it
	p->ntempmax = 20;
	p->ntemp = 1;
	p->debug = 0;
	p->pthresh = 0;
	p->nstep_adj = 100000;
	p->nprintt = 500000;
	p->keepall = 0;
	p->paranoid = 0;
	p->oldh = NULL;
	p->npre = 0;
	p->hthresh = 0.70;
	p->restart = -1;
	p->binthresh = 0;
	p->tonlyadd = 0;
	p->nkeepold = 0;
	p->sum = 0;
	p->removet = 99;
	p->tstop=0;
	p->tnorm=0;
	p->noadjust = 0;
	p->nfail1=0;
	p->nfail2=2;
	p->phthresh = 0.1;
	p->ignoreb=0;
	p->deltat=0;
    	p->gmethod=0;
	p->ttarget = 0;
	p->ttarget_harvest = 0;
	p->printpdb=1000;

	while(fgets(aux,500,fin)!=NULL)
    {
       if (!strncmp(aux,"method",6))
        if ( !sscanf(aux,"method %s",p->method) ) FatalError("Cannot read method in inputfile");

       if (!strncmp(aux,"anneal",6))
         if ( !sscanf(aux,"anneal %s",p->anneal) ) FatalError("Cannot read anneal in inputfile");

       if (!strncmp(aux,"ntmax",5))
        if ( !sscanf(aux,"ntmax %d",&p->ntempmax) ) FatalError("Cannot read ntmax in inputfile");

       if (!strncmp(aux,"nstep",5))
         if ( !sscanf(aux,"nstep %d",&p->nstep) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"preamble",8))
          if ( !sscanf(aux,"preamble %d",&p->npre) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"nsadj",5))
         if ( !sscanf(aux,"nsadj %d",&p->nstep_adj) ) FatalError("Cannot read nsadj in inputfile");

       if (!strncmp(aux,"nprint",6))
         if ( !sscanf(aux,"nprint %d",&p->nprint) ) FatalError("Cannot read nprint in inputfile");

	if (!strncmp(aux,"printpdb",8))
         if ( !sscanf(aux,"printpdb %d",&p->printpdb) ) FatalError("Cannot read printpdb in inputfile");

       if (!strncmp(aux,"nprntt",6))
          if ( !sscanf(aux,"nprntt %d",&p->nprintt) ) FatalError("Cannot read nprntt in inputfile");

       if (!strncmp(aux,"ntemp",5))
         if ( !sscanf(aux,"ntemp %d",&p->ntemp) ) FatalError("Cannot read ntemp in inputfile");

       if (!strncmp(aux,"debug",5))
          if ( !sscanf(aux,"debug %d",&p->debug) ) FatalError("Cannot read debug in inputfile");

       if (!strncmp(aux,"temperatures",12))
       {
    	   if (p->ntemp == 0) FatalError("You must specify ntemp before listing the temperatures");
    	   for (i=0;i<p->ntemp;i++)
    	   {
    		   r = fscanf(fin,"%lf %lf",&input_temp[i],&dumb);
    		   if (r==2) input_g[i] = dumb;
    		   else if (r==1) input_g[i] = 0.;
    		   else FatalError("Cannot read temperature in inputfile");
    	   }
       }

       if (!strncmp(aux,"emin",4))
         if ( !sscanf(aux,"emin %lf",&p->emin) ) FatalError("Cannot read emin in inputfile");

       if (!strncmp(aux,"emax",4))
          if ( !sscanf(aux,"emax %lf",&p->emax) ) FatalError("Cannot read emax in inputfile");

       if (!strncmp(aux,"ebin",4))
          if ( !sscanf(aux,"ebin %lf",&p->ebin) ) FatalError("Cannot read bin in inputfile");

       if (!strncmp(aux,"tfile",5))
          if ( !sscanf(aux,"tfile %s",p->nftout) ) FatalError("Cannot read tfile in inputfile");

       if (!strncmp(aux,"thefile",7))
          if ( !sscanf(aux,"thefile %s",p->nfthe) ) FatalError("Cannot read thefile in inputfile");

       if (!strncmp(aux,"dosfile",7))
          if ( !sscanf(aux,"dosfile %s",p->nfdos) ) FatalError("Cannot read dosfile in inputfile");

       if (!strncmp(aux,"histofile",9))
           if ( !sscanf(aux,"histofile %s",p->nfhisto) ) FatalError("Cannot read histofile in inputfile");

       if (!strncmp(aux,"dumbfile",8))
          if ( !sscanf(aux,"dumbfile %s",p->nfdumb) ) FatalError("Cannot read dumbfile in inputfile");

       if (!strncmp(aux,"currentdumbfile",15))
          if ( !sscanf(aux,"currentdumbfile %s",p->nfdumb2) ) FatalError("Cannot read currentdumbfile in inputfile");

       if (!strncmp(aux,"k_new",5))
           if ( !sscanf(aux,"k_new %lf",&p->k) ) FatalError("Cannot read k_new in inputfile");

       if (!strncmp(aux,"lpthresh",8))
           if ( !sscanf(aux,"lpthresh %lf",&p->pthresh) ) FatalError("Cannot read lpthresh in inputfile");

       if (!strncmp(aux,"hthresh",7))
                 if ( !sscanf(aux,"hthresh %lf",&p->hthresh) ) FatalError("Cannot read hthresh in inputfile");

       if (!strncmp(aux,"keepall",7)) p->keepall = 1;

       if (!strncmp(aux,"paranoid",8)) p->paranoid = 1;

       if (!strncmp(aux,"restart",7))
    	   if ( !sscanf(aux,"restart %d",&p->restart) ) FatalError("Cannot read restart in inputfile");

       if (!strncmp(aux,"lp_new",6))
           if ( !sscanf(aux,"lp_new %lf",&p->p_new) ) FatalError("Cannot read lp_new in inputfile");

       if (!strncmp(aux,"t_only_add",10))
                 if ( !sscanf(aux,"t_only_add %lf",&p->tonlyadd) ) FatalError("Cannot read t_only_add in inputfile");

       if (!strncmp(aux,"binthresh",9))
                        if ( !sscanf(aux,"binthresh %lf",&p->binthresh) ) FatalError("Cannot read binthresh in inputfile");

       if (!strncmp(aux,"nkeepold",8))
           if ( !sscanf(aux,"nkeepold %d",&p->nkeepold) ) FatalError("Cannot read nkeepold in inputfile");

       if (!strncmp(aux,"sumoldhisto",11))
    	   if ( !sscanf(aux,"sumoldhisto %d",&p->sum) ) FatalError("Cannot read sumoldhisto in inputfile");

       if (!strncmp(aux,"removet",7))
    	   if ( !sscanf(aux,"removet %d",&p->removet) ) FatalError("Cannot read removet in inputfile");

       if (!strncmp(aux,"tstop",5))
      	   if ( !sscanf(aux,"tstop %lf",&p->tstop) ) FatalError("Cannot read tstop in inputfile");

       if (!strncmp(aux,"ttarget",7))
      	   if ( !sscanf(aux,"ttarget %lf",&p->ttarget) ) FatalError("Cannot read ttarget in inputfile");

       if (!strncmp(aux,"tnorm",5))
      	   if ( !sscanf(aux,"tnorm %lf",&p->tnorm) ) FatalError("Cannot read tnorm in inputfile");

	if (!strncmp(aux,"nfail1",6))
      	   if ( !sscanf(aux,"nfail1 %d",&p->nfail1) ) FatalError("Cannot read nfail1 in inputfile");

	if (!strncmp(aux,"nfail2",6))
      	   if ( !sscanf(aux,"nfail2 %d",&p->nfail2) ) FatalError("Cannot read nfail2 in inputfile");

	if (!strncmp(aux,"phthresh",8))
               if ( !sscanf(aux,"phthresh %lf",&p->phthresh) ) FatalError("Cannot read phthresh in inputfile");

	if (!strncmp(aux,"ignoreb",7))
		if ( !sscanf(aux,"ignoreb %d",&p->ignoreb) ) FatalError("Cannot read ignoreb in inputfile");

	if (!strncmp(aux,"deltat",6))
   		if ( !sscanf(aux,"deltat %lf",&p->deltat) ) FatalError("Cannot read deltat in inputfile");
	if (!strncmp(aux,"gmethod",7))
   		if ( !sscanf(aux,"gmethod %s",aux_gm) ) FatalError("Cannot read gmethod in inputfile");
	if (!strncmp(aux,"pdbnf",5))
   		if ( !sscanf(aux,"pdbnf %s",p->pdbnf) ) FatalError("Cannot read pdbnf in inputfile");


    }

	if (!strcmp(p->anneal,"sigma")) p->p_new = 0;
	else if (!strcmp(p->anneal,"prob"))	p->k = 0;
	

	// restart
	if (p->restart > -1)
		{
			Restart(p);
			if (p->debug>0) fprintf(stderr,"Start the simulation...\n\n");
			return p;
		}

	if (p->debug>0)
	{
		fprintf(stderr,"\n************************************\n");
		fprintf(stderr,"* SIMULATED TEMPERING              *\n");
		fprintf(stderr,"************************************\n\n");
		fflush(stderr);

		if (!strcmp(p->method,"stempering") || !strcmp(p->method,"adaptive"))
			fprintf(stderr,"method =\t%s\n",p->method);
		else FatalError("Method not known");
		fprintf(stderr,"ntempmax = \t%d\n",p->ntempmax);
		fprintf(stderr,"ntemp = \t%d\n",p->ntemp);
		fprintf(stderr,"temperatures:\n");
		for (i=0;i<p->ntemp;i++)
			fprintf(stderr,"T[%d] = %lf\n",i,input_temp[i]);
		fprintf(stderr,"weights:\n");
		for (i=0;i<p->ntemp;i++)
			fprintf(stderr,"g[%d] = %lf\n",i,input_g[i]);
		fprintf(stderr,"nstep = \t%d\n",p->nstep);
		if (!strcmp(p->method,"adaptive"))
				fprintf(stderr,"nsadj = \t%d\n",p->nstep_adj);
		if (p->npre>0) fprintf(stderr,"preamble = \t%d\n",p->npre);
		fprintf(stderr,"nprint = \t%d\n",p->nprint);
		fprintf(stderr,"nprintt = \t%d\n",p->nprintt);
		fprintf(stderr,"debug = \t%d\n",p->debug);
		fprintf(stderr,"emin = \t\t%lf\n",p->emin);
		fprintf(stderr,"emax = \t\t%lf\n",p->emax);
		fprintf(stderr,"ebin = \t\t%lf\n",p->ebin);
		fprintf(stderr,"nbin = \t\t%d\n",(int)((p->emax-p->emin)/p->ebin));
		fprintf(stderr,"binthresh = \t%lf\n",p->binthresh);
		fprintf(stderr,"phthresh = \t%lf\n",p->phthresh);
		if (p->tstop>0) fprintf(stderr,"tstop = \t\t%lf\n",p->tstop);
		if (p->tnorm>0) fprintf(stderr,"tnorm = \t\t%lf\n",p->tnorm);

		if (!strcmp(p->anneal,"sigma"))
		{
			fprintf(stderr,"Anneal = \tproportionally to energy stdev\n");
			fprintf(stderr,"k_new = \t%lf\n",p->k);
		}
		else if (!strcmp(p->anneal,"prob"))
		{
			fprintf(stderr,"Anneal = \tsetting wished exchange probability\n");
			fprintf(stderr,"p_new = \t%lf\n",p->p_new);
		}
		else
			fprintf(stderr,"Anneal = \tnone (fixed range of temperatures\n");

		if (p->removet>0)
		{
			if (p->pthresh<0)
				fprintf(stderr,"lpthresh = \t%lf (to reduce the number of temperatures)\n",p->pthresh);
			else if (p->pthresh>8)
				fprintf(stderr,"lpthresh = automatically calculated\n");
			fprintf(stderr,"removet = \t%d\n",p->removet);
		}
		fprintf(stderr,"hthresh = \t%lf (to keep a histogram)\n",p->hthresh);

		if (p->keepall==1 && p->nkeepold==0) fprintf(stderr,"Calculate g(E) from all past history\n");
		else if (p->keepall==1) fprintf(stderr,"Calculate g(E) from last %d histograms\n",p->nkeepold);
		if (p->keepall==1 && p->sum==1) fprintf(stderr,"Sum together histograms with same temperature\n");
		if (p->keepall==1 && p->sum==2) fprintf(stderr,"Delete histograms with same temperature than newer\n");

		if (p->paranoid==1) fprintf(stderr,"Be paranoid\n");
		if (p->tonlyadd>0) fprintf(stderr,"Below T=%lf, only add temperatures\n",p->tonlyadd);
		if (strcmp(p->nfhisto,"")) fprintf(stderr,"Write histograms to %s\n",p->nfhisto);
	}

	// Check the allocation dimensions
	if (NHISTOMAX>NTEMPMAX) FatalError("NHISTOMAX should be smaller than NTEMPMAX");
	if (p->ntempmax>NTEMPMAX) FatalError("ntmax should be smaller than NTEMPMAX");
	if (p->ntemp>p->ntempmax) FatalError("ntemp should be smaller than ntmax");

	// Initialize common stuff
	p->ftout = fopen(p->nftout,"w");
	if (p->ftout==NULL) FatalError("Error opening 'tfile' for writing");
	fprintf(p->ftout,"# STEP\tTEMP\tNTEMP\n");
	p->f = NULL;
	p->nbin = (int) ((p->emax-p->emin)/p->ebin) + 1;
	p->h = AlloDoubleMat(p->ntempmax,p->nbin);
	p->binok = AlloInt(p->nbin);
	p->htmp = AlloDoubleMat(NHISTOMAX,p->nbin);
    fprintf(stderr,"allocate htmp %dx%d\n",NHISTOMAX,p->nbin);

	p->out = AlloDoubleMat(4,NTBIN);
	p->boltzp = AlloDoubleMat(NHISTOMAX,p->nbin);

	// Initialize standard simulated tempering
	if (!strcmp(p->method,"stempering"))
	{
		if (p->ntemp<2) FatalError("You need more than two tempering to make a stempering");
		p->nm = 1;
		p->temp = AlloDouble(p->ntemp);
		p->g = AlloDouble(p->ntemp);
		p->prob_down = AlloDouble(p->ntemp);
		p->prob_up = AlloDouble(p->ntemp);
		p->counts = AlloInt(p->ntemp);

		for (i=0;i<p->ntemp;i++)
		{
			p->temp[i] = input_temp[i];
			p->g[i] = input_g[i];
		}
	}

	// Initialize adaptive simulated tempering
	else if (!strcmp(p->method,"adaptive"))
	{
		p->nm = 2;
		p->temp = AlloDouble(p->ntempmax);
		p->g = AlloDouble(p->ntempmax);
		p->prob_down = AlloDouble(p->ntempmax);
		p->prob_up = AlloDouble(p->ntempmax);
		p->counts = AlloInt(p->ntempmax);
		p->f = AlloDouble(p->ntempmax);
		p->reliable_lg = AlloDouble(p->nbin);
		p->current_lg = AlloDouble(p->nbin);
		p->reliable_t = AlloDouble(p->ntempmax);
		p->reliable_g = AlloDouble(p->ntempmax);

		for (i=0;i<p->ntemp;i++)
			p->g[i] = input_g[i];

		if (p->keepall==1)
		{
			p->oldh = AlloDoubleMat(NHISTOMAX,p->nbin);
			p->oldt = AlloDouble(NHISTOMAX);
			p->oldt_iter = AlloInt(NHISTOMAX);
			p->noldt = 0;

			if (p->debug>0) fprintf(stderr,"Allocate memory to keep all histograms up to %d.\n",NHISTOMAX);
		}

		p->fpacc = fopen("accept.dat","w");
		if (!p->fpacc) FatalError("cannot open accetp.dat");

	}

	for (i=0;i<p->ntemp;i++)
			p->temp[i] = input_temp[i];

	// structure for restart
	p->st_restart = AlloRestart(NTEMPMAX,p->nbin,NRESTARTS);

	// Start from highest temperature
	p->itemp = 0;

	// Reset counters
	p->count_st = 0;
	p->count_adj = 0;
	p->count_print = 0;
	p->count_printt = 0;
	p->iter = 0;
	p->failure = 0;

	return p;
}

/*************************************************************
 * At each step of the simulation, before moving, collect
 * histograms, decide whether to attempt temperature changes
 * and do it.
 * 		it needs the current potential energy
 * 		the new temperature is in (double) p->temp[p->itemp]
 * 		returns 1 if change is accepted, 0 if not.
 *************************************************************/
int STempering(double energy, double step, struct st_stuff *p)
{
	int ie,r=0;
	double tmin, tbin, curr_temp;
	double *lg=NULL;
	FILE *fthe;

	// collect histograms
	if (energy<p->emin || energy>p->emax)		// if energy out of range
	{
			fprintf(stderr,"\n\nE=%lf (emin=%lf emax=%lf)\n",energy,p->emin,p->emax);
			if (p->paranoid==0) fprintf(stderr,"WARNING: energy out of range\n");
			else FatalError("Energy out of range");
	}
 	else if (p->npre==0 || step>(double)p->npre)
    	{
            // if energy within range
		ie = (int)( (double) (energy - p->emin) / p->ebin + 0.5); 	//NB: since (energy - p->emin) / p->ebin > 0 to round properly I just have to add 0.5
			p->h[p->itemp][ie] ++;
	}


	/// make simulated tempering
	if (p->count_st >= p->nstep && step > 0 && p->ntemp > 1)				// every nstep
		if (p->npre==0 || step > (double)p->npre)					// skip first steps
		{
			r = DoSTempering(energy,p);
			p->count_st=0;
		}

	// print temperature
    	if (p->count_print >= p->nprint)
   	{
    		fprintf(p->ftout,"%.0lf\t%lf\t%d\n",step,p->temp[p->itemp],p->ntemp);
    		fflush(p->ftout);
    		p->count_print = 0;
    	}


	// if not adaptive, print thermodynamics, dumb energies, statistics and histograms
	if (p->nm==1 && p->count_printt>(double)p->nprintt)
	{
		// histograms
		if (strcmp(p->nfhisto,""))
				PrintHistogram(p->nfhisto,p->h,p->temp,p->ntemp,p->nbin,p->emin,p->ebin);

		// statistics
		PrintStatistics(p, step);

		// dumb energies
		fthe = fopen(p->nfdumb,"w");
		if (!fthe) FatalError("Cannot open file for writing dumb averages");
		PrintAverageEnergies(fthe,p->h,p->temp,p->ntemp,p->emin,p->ebin,p->nbin);
		fclose(fthe);

		// thermodynamics
		if (strcmp(p->nfthe,""))
		{
			if (p->debug>0) fprintf(stderr,"Calculates thermdynamics\n");
			OrderTemperatures(p->ntemp,p->temp,p->g,p->prob_up,p->prob_down,p->h,p->nbin,0,p->counts);
			FilterHistograms(p->h,p->htmp,p->ntemp,p->nbin,p->binthresh);
			lg = MHistogram(p->htmp,p->temp,p->ntemp,p->nbin,p->emin,p->ebin,0,KB,lg,p->f,0,p->ignoreb,NULL,p->deltat,p->binok);
			fthe = fopen(p->nfthe,"w");
			if (!fthe) FatalError("Cannot open file for writing thermodynamics");
			tmin = p->temp[p->ntemp-1]/2;
			tbin = (p->temp[0]*3/2 - tmin)/NTBIN;
			CalculateThermodynamics(lg,p->emin,p->ebin,p->nbin,tmin,tbin,NTBIN,p->out,KB,NULL,NULL,p->binok);
			PrintThermodynamics(fthe,p->out,tmin,tbin,NTBIN);
			fclose(fthe);
		}

		p->count_printt = 0;
	}

	// adjust temperatures
	if (p->npre==0 || step>(double)p->npre)						// skip first steps
		if (p->nm == 2)
			if ((p->count_adj >= p->nstep_adj*p->ntemp ) && step>0) 	// every nstep_adj*ntemp
			{
				curr_temp = p->temp[p->itemp];				//save current T

				AdjustSTempering(p,step);

				p->itemp = FindClosestT(curr_temp,p->ntemp,p->temp);	//restart from closest T
											//(to avoid problems if the current T has been deleted for instance)
				r = 1;
				p->count_adj = 0;					// reset all counters


				// below this temperature, stop
				if (p->tnorm>0)
					if (p->temp[p->ntemp-1]<p->tnorm)
					{
						if (p->debug>0) fprintf(stderr,"Reached target temperature (Tstop=%lf).\nCONTINUING WITH STEMPERING\n",p->tnorm);
                        			p->nm=1;
					}
				// below this temperature, set to ttarget and make a last iteration
				
				if (p->ttarget>0 )
					if (p->temp[p->ntemp-1]<=p->ttarget)
					{
					//MODIFICHE: esco solo quando ho stampato il numero di conformazioni stabilite
					////////if (p->ttarget_harvest)
					////////{
					////////	if (p->debug>0) fprintf(stderr,"Reached target temperature (Ttarget=%lf).\nENDING SIMULATION\n",p->ttarget);
					////////	exit(0);
					////////}

						p->temp[p->ntemp-1]=p->ttarget;
						p->g[p->ntemp-1] = CalculateG(p->current_lg, p->emin, p->ebin, p->nbin, p->temp[p->ntemp-1], KB, p->binok, p->gmethod) - 
								CalculateG(p->current_lg, p->emin, p->ebin, p->nbin, p->temp[0], KB, p->binok, p->gmethod);
						if (p->debug>0) fprintf(stderr,"Reached target temperature (Ttarget=%lf).\nHARVESTING RUN\n",p->ttarget);
									
						p->ttarget_harvest = 1;	
						p->pdbf = fopen(p->pdbnf,"w");
						if (!p->pdbf) FatalError("Cannot open pdbnf file for harvesting");
					}

			}

	if (p->tstop>0)
		if (p->temp[p->ntemp-1]<p->tstop)
		{
			if (p->debug>0) fprintf(stderr,"Reached target temperature (Tstop=%lf).\nENDING SIMULATION\n",p->tstop);
			exit(0);
		}

	// advance clocks
	p->count_st ++;
	p->count_adj ++;
	p->count_print ++;
	//MODIFICHE: mi serve solo se ho static simuled tempering
	if(p->nm==1)
	p->count_printt ++;

        if (step<2) r=1;		// at the first step, remember to set the temperatures defined in _STEMPERING

	return r;
}

/*************************************************************
 * Change temperature
 * 		returns 1 if change is accepted, 0 if not.
 *************************************************************/
int DoSTempering(double energy, struct st_stuff *p)
{
	double beta1,beta2,delta;
	int iw,acc;

	if (p->debug>1) fprintf(stderr,"\nStempering: change temperature\n");
	if (p->debug>2) fprintf(stderr,"Current T=%lf (itemp=%d)\n",p->temp[p->itemp],p->itemp);

	// update statistics
    p->counts[p->itemp] ++;

    // decide whether to go up or down
    if (p->itemp==0) iw=1;
    else if (p->itemp==p->ntemp-1) iw=-1;
    else
    {
    	if (rand()>RAND_MAX/2) iw=1;
    	else iw=-1;
    }

	// evaluate jump probablity
    beta1 = 1./(KB*p->temp[p->itemp]);
    beta2 = 1./(KB*p->temp[p->itemp+iw]);
    delta = (beta1-beta2) * energy + p->g[p->itemp+iw] - p->g[p->itemp];


    if (delta>0 || (double)rand()/RAND_MAX<exp(delta))
    {
		if (iw == 1) p->prob_down[p->itemp]++;
		else p->prob_up[p->itemp]++;
		p->itemp += iw;
		acc = 1;
    }
    else acc = 0;


	if (p->debug>1)
	{
		if (acc==1) fprintf(stderr,"new T=%lf\n",p->temp[p->itemp]);
		else fprintf(stderr,"mantaining old temperature (T=%lf)\n",p->temp[p->itemp]);
		fflush(stderr);
	}

	return acc;
}

/*****************************************************************************
 Reorder temperatures starting from highest
 *****************************************************************************/
void OrderTemperatures(int ntemp, double *temp, double *g, double *prob_up, double *prob_down, double **h,int nbin, int debug, int *counts)
{
  int i,j,k;
  double d;

  for (i=0;i<ntemp-1;i++)
    for (j=i+1;j<ntemp;j++)
      if ( temp[i] < temp[j] )
	  {
		 d = temp[i];
		 temp[i] = temp[j];
		 temp[j] = d;

		 if (g != NULL)
		 {
			 d = g[i];
			 g[i] = g[j];
			 g[j] = d;
		 }

		 if (prob_up != NULL)
		 {
			 d = prob_up[i];
			 prob_up[i] = prob_up[j];
			 prob_up[j] = d;
		 }

		 if (prob_down != NULL)
		 {
			 d = prob_down[i];
			 prob_down[i] = prob_down[j];
			 prob_down[j] = d;
		 }

		 if (h != NULL)
			 for (k=0;k<nbin;k++)
			 {
				 d = h[i][k];
				 h[i][k] = h[j][k];
				 h[j][k] = d;
			 }

		 if (counts != NULL)
		 {
			 k = counts[i];
			 counts[i] = counts[j];
			 counts[j] = k;
		 }
	  }

  if (debug>2) {
	fprintf(stderr,"\nReordering temperatures\n new order: ");
	for (i=0;i<ntemp;i++) fprintf(stderr,"%lf ",temp[i]);
	fprintf(stderr,"\n");
 }
}

/*****************************************************************************
 Calculates the average energy from a straightforward average of the histogram
 *****************************************************************************/
void PrintAverageEnergies(FILE *fout, double **h, double *temp, int ntemp, double emin, double ebin, int nbin)
{
	int it,ie;
	double em,z,e,s;

	fprintf(fout,"#TEMP\t<E>\tsigma_E\tNormalization\n");
	for (it=0;it<ntemp;it++)
	{
		em = z = s = 0.;

		for (ie=0;ie<nbin;ie++)
		{
			e = (double) ie * ebin + emin;
			em += e * h[it][ie];
			s += (e * e * h[it][ie]);
			z += h[it][ie];
		}
		em = em/z;
		fprintf(fout,"%lf\t%lf\t%lf\t%lf\n",temp[it],em,sqrt(s/z-em*em),z);
	}
}

void PrintHistogram(char *filename, double **h, double *t, int ntemp, int nbin, double emin, double ebin)
{
	int ie,it;
	double e;
	FILE *fh;


	fh = fopen(filename,"w");
	if (!fh) FatalError("Cannot open histogram file");
	fprintf(fh,"%%E\t\t");
	for (it=0;it<ntemp;it++) fprintf(fh,"T=%lf\t",t[it]);
	fprintf(fh,"\n");
	for (ie=0;ie<nbin;ie++)
	{
		e = (double) ie * ebin + emin;
		fprintf(fh,"%.8lf\t",e);
		for (it=0;it<ntemp;it++) fprintf(fh,"%.6e\t",h[it][ie]);
		fprintf(fh,"\n");
	}
	fclose(fh);

}

void Restart(struct st_stuff *p)
{
	int i,daux,w=-1;
	double faux1,faux2,step;
	FILE *fr;
	char aux[100];

	fr = fopen("RESTART_ST","r");
	if (!fr) FatalError("Cannot open restart file");

	if (p->debug>0)
	{
		fprintf(stderr,"\n************************************\n");
		fprintf(stderr,"* SIMULATED TEMPERING              *\n");
		fprintf(stderr,"************************************\n\n");
		fprintf(stderr,"\nRESTART from RESTART_ST (snap=%d)\n\n",p->restart);
		fflush(stderr);
	}

	// find the right snapshot to restart
	do
	{
		w = -1;
		if (fgets(aux,100,fr)==NULL) break;
		sscanf(aux,"restart %d",&w);
	} while (w != p->restart);

	if (w==-1) FatalError("Cannot find the requested snapshot in RESTART");

	// read ntemp
	if (fscanf(fr,"%lf %d",&step,&p->ntemp)!=2) FatalError("Wrong format in restart file [1]");

	fprintf(stderr,"Restarting from step %lf of the preceding simulation.\n",step);

	// allocate stuff for standard stempering
	if (!strcmp(p->method,"stempering"))
	{
		if (p->ntemp<2) FatalError("You need more than two tempering to make a stempering");
		p->nm = 1;
		p->temp = AlloDouble(p->ntemp);
		p->g = AlloDouble(p->ntemp);
		p->prob_down = AlloDouble(p->ntemp);
		p->prob_up = AlloDouble(p->ntemp);
		p->counts = AlloInt(p->ntemp);
		p->keepall = 0;
	}

	// allocate staff for continuing adaptive stempering
	else if (!strcmp(p->method,"adaptive"))
	{
		p->nm = 2;
		p->temp = AlloDouble(p->ntempmax);
		p->g = AlloDouble(p->ntempmax);
		p->prob_down = AlloDouble(p->ntempmax);
		p->prob_up = AlloDouble(p->ntempmax);
		p->counts = AlloInt(p->ntempmax);
		p->f = AlloDouble(p->ntempmax);
		p->reliable_t = AlloDouble(p->ntempmax);
		p->reliable_g = AlloDouble(p->ntempmax);

		p->fpacc = fopen("accept.dat","w");
		if (!p->fpacc) FatalError("cannot open accept.dat");

	}

	// read temperatures, weights
	for (i=0;i<p->ntemp;i++) if(fscanf(fr,"%lf %lf",&p->temp[i],&p->g[i])!=2) FatalError("Wrong format in restart file [2]");
	if (p->debug>0)
	{
		fprintf(stderr,"Read %d temperatures:\n",p->ntemp);
		for (i=0;i<p->ntemp;i++) fprintf(stderr,"T=%lf\tg=%lf\n",p->temp[i],p->g[i]);
	}

	// if you want to change them, read binning (override STEMPERING.dat)
	if(fscanf(fr,"%lf %lf %d",&faux1,&faux2,&daux)==3)
	{
		p->emin = faux1;
		p->ebin = faux2;
		p->nbin = daux;
		p->emax = p->emin + p->ebin * p->nbin;
		if (p->debug>0) fprintf(stderr,"Binning read from restart\n");
	}
	else if (p->debug>0) fprintf(stderr,"Binning read from STEMPERING\n");

	// allocate density of states
	p->reliable_lg = AlloDouble(p->nbin);
	p->current_lg = AlloDouble(p->nbin);

	// if you want, read density of states
	if(fscanf(fr,"%lf ",&p->reliable_lg[0])==1)
	{
		for (i=1;i<p->nbin;i++)
			if(fscanf(fr,"%lf ",&p->reliable_lg[i])!=1) FatalError("Wrong format in restart file [5]");
		if (p->debug>0) fprintf(stderr,"Density of states read from restart file\n");
	}


	if (p->debug>0) fprintf(stderr,"Restart file read correctly!\n");
	fclose(fr);


	if (p->debug>0)
	{
		if (!strcmp(p->method,"stempering") || !strcmp(p->method,"adaptive"))
			fprintf(stderr,"method =\t%s\n",p->method);
		else FatalError("Method not known");
		fprintf(stderr,"ntempmax = \t%d\n",p->ntempmax);
		fprintf(stderr,"nstep = \t%d\n",p->nstep);
		if (!strcmp(p->method,"adaptive"))
				fprintf(stderr,"nsadj = \t%d\n",p->nstep_adj);
		if (p->npre>0) fprintf(stderr,"preamble = \t%d\n",p->npre);
		fprintf(stderr,"nprint = \t%d\n",p->nprint);
		fprintf(stderr,"nprintt = \t%d\n",p->nprintt);
		fprintf(stderr,"debug = \t%d\n",p->debug);
		fprintf(stderr,"emin = \t\t%lf\n",p->emin);
		fprintf(stderr,"emax = \t\t%lf\n",p->emax);
		fprintf(stderr,"ebin = \t\t%lf\n",p->ebin);
		fprintf(stderr,"nbin = \t\t%d\n",(int)((p->emax-p->emin)/p->ebin));
		fprintf(stderr,"binthresh = \t%lf\n",p->binthresh);

		if (!strcmp(p->anneal,"sigma"))
		{
			fprintf(stderr,"Anneal = \tproportionally to energy stdev\n");
			fprintf(stderr,"k_new = \t%lf\n",p->k);
			p->p_new = 0;
		}
		else if (!strcmp(p->anneal,"prob"))
		{
			fprintf(stderr,"Anneal = \tsetting wished exchange probability\n");
			fprintf(stderr,"p_new = \t%lf\n",p->p_new);
			p->k = 0;
		}
		else
			fprintf(stderr,"Anneal = \tnone (fixed range of temperatures\n");

		if (p->removet>0)
		{
			if (p->pthresh<0)
				fprintf(stderr,"lpthresh = \t%lf (to reduce the number of temperatures)\n",p->pthresh);
			else if (p->pthresh>8)
				fprintf(stderr,"lpthresh = automatically calculated\n");
			fprintf(stderr,"removet = \t%d\n",p->removet);
		}
		fprintf(stderr,"hthresh = \t%lf (to keep a histogram)\n",p->hthresh);

		if (p->keepall==1 && p->nkeepold==0) fprintf(stderr,"Calculate g(E) from all past history\n");
		else if (p->keepall==1) fprintf(stderr,"Calculate g(E) from last %d histograms\n",p->nkeepold);
		if (p->keepall==1 && p->sum==1) fprintf(stderr,"Sum together histograms with same temperature\n");
		if (p->keepall==1 && p->sum==2) fprintf(stderr,"Delete histograms with same temperature than newer\n");

		if (p->paranoid==1) fprintf(stderr,"Be paranoid\n");
		if (p->tonlyadd>0) fprintf(stderr,"Below T=%lf, only add temperatures\n",p->tonlyadd);
		if (strcmp(p->nfhisto,"")) fprintf(stderr,"Write histograms to %s\n",p->nfhisto);
	}

	// Check the allocation dimensions
	if (NHISTOMAX>NTEMPMAX) FatalError("NHISTOMAX should be smaller than NTEMPMAX");
	if (p->ntempmax>NTEMPMAX) FatalError("ntmax should be smaller than NTEMPMAX");
	if (p->ntemp>p->ntempmax) FatalError("ntemp should be smaller than ntmax");

	// Initialize common stuff
	p->ftout = fopen(p->nftout,"w");
	if (p->ftout==NULL) FatalError("Error opening 'tfile' for writing");
	fprintf(p->ftout,"# STEP\tTEMP\tNTEMP\n");
	if (p->debug>0) fprintf(stderr,"Open %s for writing temperatures\n",p->nftout);
	p->binok = AlloInt(p->nbin);

	p->f = NULL;

	// Initialize other stuff
	p->h = AlloDoubleMat(p->ntempmax,p->nbin);
	p->htmp = AlloDoubleMat(NHISTOMAX,p->nbin);
	fprintf(stderr,"allocate htmp %dx%d\n",NHISTOMAX,p->nbin);
	p->out = AlloDoubleMat(4,NTBIN);
	p->boltzp = AlloDoubleMat(NHISTOMAX,p->nbin);

	if (p->keepall==1)
	{
		p->oldh = AlloDoubleMat(NHISTOMAX,p->nbin);
		p->oldt = AlloDouble(NHISTOMAX);
		p->oldt_iter = AlloInt(NHISTOMAX);
		p->noldt = 0;

		if (p->debug>0) fprintf(stderr,"Allocate memory to keep all histograms up to %d.\n",NHISTOMAX);
	}
	else p->oldh = NULL;

	// structure for restart
	p->st_restart = AlloRestart(NTEMPMAX,p->nbin,NRESTARTS);

	// Start from highest temperature
	p->itemp = 0;

	// Reset counters
	p->count_st = 0;
	p->count_adj = 0;
	p->count_print = 0;
	p->count_printt = 0;
	p->iter = 0;
	p->failure = 0;
}

int FindClosestT(double t, int ntemp, double *temp)
{
	int i,itemp=0;
	double mindist=9999.;

	for (i=0;i<ntemp;i++)
	{
		if (fabs(temp[i]-t)<mindist)
		{
			mindist = fabs(temp[i]-t);
			itemp = i;
		}
	}

	if (mindist>9998.) FatalError("Cannot find closest temperature");

	return itemp;
}

void PrintStatistics(struct st_stuff *p, double step)
{
	int i;
	double p_down,p_up,p_tot;

	fprintf(stderr,"\n\nSIMULATED TEMPERING STATISTICS  (step=%.0lf):\n",step);
	p_tot = 1.0;
	for (i=0;i<p->ntemp;i++)
	{
			p_down = p_up = 0.;
			if (i==0) {p_down = p->prob_down[i]/p->counts[i]; p_tot *= p_down;}
			else if (i==p->ntemp-1) {p_up = p->prob_up[i]/p->counts[i]; p_tot *= p_up;}
			else
			{
				p_up = p->prob_up[i]/p->counts[i];
				p_down = p->prob_down[i]/p->counts[i];
				p_tot *= p_up * p_down;
    		}
			fprintf(stderr,"\tT=%lf\tp_up=%lf\tp_down=%lf\tg=%lf\n",p->temp[i],p_up,p_down,p->g[i]);
	}
	fprintf(stderr,"lp_tot/n=%lf\n",log(p_tot)/p->ntemp);
	fflush(stderr);
}
