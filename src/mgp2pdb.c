/*
 * mgp2pdb.c
 *
 *  Created on: Jun 26, 2016
 *      Author: riccardo
 */

#include "montegrappa.h"
#include <ctype.h>
#include <unistd.h>

void mgp2pdbHeader();
int polymercounter(char* filename);

int main(int argc, char *argv[])
{
	struct s_polymer *polymer;
	char pol_file[100] = " ", pdb_file[100] = " ";
	
	int npol=1;
	int i,nat,ntypes;
	
	FILE *logfile=NULL, *pdbout=NULL;
	
	mgp2pdbHeader();
	
	logfile = fopen("mgp2pdb.log","w");
	if (argc == 1)
	{
		fprintf(stderr,"usage: mgp2pdb -i polfile -o pdbfile\n");
		return 1;
	}
	
	
	
	// Options handling
	while ((i = getopt (argc, argv, "hi:o:")) != -1)
		switch (i)
		{
			case 'h':
				fprintf(stderr,"usage: mgp2pdb -i polfile -o pdbfile\n");
				return 1;
				break;
			case 'i':
				strcpy(pol_file,optarg);
				fprintf(stderr," -> input pol file: %s\n",pol_file);
				fprintf(logfile,"Input pol file: %s\n",pol_file);
				break;
			case 'o':
				strcpy(pdb_file,optarg);
				fprintf(stderr," -> output pdb file: %s\n",pdb_file);
				fprintf(logfile,"Output pdb file: %s\n",pdb_file);
				break;
			case '?':
				if (optopt == 'i')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (optopt == 'o')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr,
							 "Unknown option character `\\x%x'.\n",
							 optopt);
				return 1;
			default:
				abort ();
	}
	
	npol = polymercounter(pol_file);
	
	fprintf(stderr," -> number of chains: %d\n",npol);
	fprintf(logfile,"Number of chains: %d\n",npol);

	
	pdbout = fopen(pdb_file,"w");
	
	
	// Polymer allocation
	polymer = AlloPolymer(npol,NRESMAX,NSIDEMAX,NROTMAX,NATOMMAX,0,0,logfile);
	
	// Polymer reading
	ReadPolymer(pol_file,polymer,logfile,npol,0, &nat, &ntypes);

	// SideChains adding
	fprintf(stderr," -> adding sidechains\n");
	for (i=0;i<npol;++i)
		AddSidechain(polymer,0,polymer->nback,i);
	
	// Polymer printing
	fprintf(logfile,"Printing pdb file %s\n",pdb_file);
	fprintf(stderr," -> printing %s\n",pdb_file);
	PrintPDBStream(polymer,npol,pdbout);
	
	
	FreePolymer(polymer,npol,NRESMAX,NSIDEMAX,0,0);

	fclose(pdbout);
	fclose(logfile);
	
	
	
	return 0;
}

void mgp2pdbHeader(){
	fprintf(stderr,"\n\n");
	fprintf(stderr,"                       ___            _ _ \n");
	fprintf(stderr,"                      |__ \\          | | |    \n");
	fprintf(stderr,"  _ __ ___   __ _ _ __   ) |_ __   __| | |__  \n");
	fprintf(stderr," | '_ ` _ \\ / _` | '_ \\ / /| '_ \\ / _` | '_ \\ \n");
	fprintf(stderr," | | | | | | (_| | |_) / /_| |_) | (_| | |_) |\n");
	fprintf(stderr," |_| |_| |_|\\__, | .__/____| .__/ \\__,_|_.__/\n");
	fprintf(stderr,"             __/ | |       | |\n");
	fprintf(stderr,"            |___/|_|       |_|\n\n\n");
	fprintf(stderr,"\t R. Capelli, 2016\n\n");
}

int polymercounter(char* filename){
	int npol = 0, chain=-1, old_chain=-1, dummy_d;
	char dummy_s[10], aux[500];
	double dummy_lf;
	
	FILE *fp = NULL;
	
	fp = fopen(filename,"r");
	while(fgets(aux,500,fp)!=NULL){
		if(sscanf( aux, "%d %d %s %d %s %d %d %lf %lf %lf %d", \
				  &dummy_d,&dummy_d,dummy_s,&dummy_d,dummy_s,&dummy_d,&chain,&dummy_lf,&dummy_lf,&dummy_lf,&dummy_d) == 11 )
		{
			if (old_chain == chain)
			{
				continue;
			}
			else
			{
				++npol;
			}
			old_chain = chain;
		}
	}
	fclose(fp);
	

	return npol;
}

